//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
/// \file DetectorConstruction.cc
/// \brief DetectorConstruction 类的实现：构建细胞微剂量学几何
///
/// 细胞几何为同心球嵌套结构：World(水盒) -> Cell(细胞膜边界) -> Nucleus(细胞核)。
/// 敏感探测器(SD)置于细胞核，域(site)采样与核级打分均在核内进行。
/// 当前所有体积材料均为液态水（细胞核成分后续可进一步细化）。

#include "DetectorConstruction.hh"

#include "G4Box.hh"
#include "G4Orb.hh"
#include "G4Exception.hh"
#include "G4Material.hh"
#include "G4NistManager.hh"
#include "G4PVPlacement.hh"
#include "G4Region.hh"
#include "G4RunManager.hh"
#include "G4SDManager.hh"
#include "G4StateManager.hh"
#include "G4SystemOfUnits.hh"
#include "globals.hh"

#include "DetectorMessenger.hh"
#include "TrackerSD.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::DetectorConstruction() : G4VUserDetectorConstruction()
{
  // 创建探测器交互命令对象，用于在运行时通过 UI 命令修改探测器参数
  fDetectorMessenger = std::make_unique<DetectorMessenger>(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::~DetectorConstruction() = default;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::Construct()
{
  /// Geant4 入口：构建整个探测器几何。
  /// 顺序为「定义材料 -> 定义体积」，体积采用同心球嵌套结构。
  /// @return 世界物理体积指针（几何树根节点）

  // 第一步：定义所需的介质材料（当前为水）
  DefineMaterials();

  // 第二步：定义体积层次（世界内含细胞，细胞内含细胞核）
  G4VPhysicalVolume* world = DefineWorld();

  return world;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::DefineMaterials()
{
  /// 定义模拟所用的介质材料。
  /// 当前直接从 NIST 数据库取用液态水（G4_WATER），后续可在此细化核成分。

  // 获取 NIST 材料管理器单例
  G4NistManager* nist = G4NistManager::Instance();

  // 若尚未指定材料，则默认使用 NIST 液态水
  if (!fMat) {
    fMat = nist->FindOrBuildMaterial("G4_WATER");
  }

  // 打印当前材料名称
  G4cout << "材料: " << fMat->GetName() << G4endl;
  // 打印完整的材料表（Geant4 内置输出，便于核对材料定义）
  G4cout << *(G4Material::GetMaterialTable()) << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::CheckConsistency()
{
  /// 校验几何参数的一致性与合法性。
  /// 任一条件不满足，都会抛出致命异常(FatalException)终止运行。

  // 检查 1：次级电子最大射程必须为正
  if (!(fMaxRange > 0.)) {
    G4ExceptionDescription msg;
    msg << "fMaxRange 必须 > 0。\n"
        << "建议参考 Tabata 公式估算次级电子的最大射程：\n"
        << "\t https://doi.org/10.1016/0029-554X(72)90463-6" << G4endl
        << "说明：此参数决定世界体积在细胞外的水层厚度（用于维持次级电子平衡）。"
        << G4endl;
    G4Exception("DetectorConstruction::CheckConsistency()", "DetCons0001",
                FatalException, msg);
  }

  // 检查 2：细胞核半径必须为正
  if (!(fNucleusRadius > 0.)) {
    G4ExceptionDescription msg;
    msg << "fNucleusRadius 必须 > 0。\n"
        << "细胞核半径必须为正，否则无法在核内放置采样位点。" << G4endl;
    G4Exception("DetectorConstruction::CheckConsistency()", "DetCons0002",
                FatalException, msg);
  }

  // 检查 3：细胞半径必须大于细胞核半径（同心球结构要求）
  if (!(fCellRadius > fNucleusRadius)) {
    G4ExceptionDescription msg;
    msg << "fCellRadius 必须 > fNucleusRadius。\n"
        << "细胞半径必须大于细胞核半径（同心球结构）。" << G4endl;
    G4Exception("DetectorConstruction::CheckConsistency()", "DetCons0003",
                FatalException, msg);
  }

  // 检查 4：域(site)半径必须为正且不大于细胞核半径
  if (!(fSiteRadius > 0. && fSiteRadius <= fNucleusRadius)) {
    G4ExceptionDescription msg;
    msg << "fSiteRadius 必须 > 0 且 <= fNucleusRadius。\n"
        << "域(site)半径必须为正且不大于细胞核半径。" << G4endl;
    G4Exception("DetectorConstruction::CheckConsistency()", "DetCons0004",
                FatalException, msg);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::DefineWorld()
{
  /// 构建世界体积（最外层水盒），并在其内部依次构建细胞与细胞核。
  /// @return 世界物理体积指针，作为整个几何树的根节点

  // 构建前先校验几何参数一致性，避免出现非法尺寸导致穿模
  CheckConsistency();

  // 世界为正方体水盒：半边长 = 细胞半径 + 次级电子最大射程
  // 这样细胞外保留足够厚的水层，用于维持次级电子平衡
  G4double worldHalf = fCellRadius + fMaxRange;  // 世界半边长

  // —— 世界几何体（正方体实体）——
  auto* solidWorld = new G4Box("World", worldHalf, worldHalf, worldHalf);

  // —— 世界逻辑体积（绑定材料，无磁场、无敏感标记）——
  auto* logicalWorld =
    new G4LogicalVolume(solidWorld, fMat, "World", nullptr, nullptr, nullptr);

  // —— 世界物理体积（置于原点，无旋转；作为根节点，母体积为空）——
  G4VPhysicalVolume* world = new G4PVPlacement(
    nullptr,            // 母体积旋转矩阵
    G4ThreeVector(),    // 相对母体积的平移（原点）
    logicalWorld,       // 关联的逻辑体积
    "World",            // 体积名称
    nullptr,            // 母逻辑体积（世界为根，故为空）
    false,              // 是否复制体积（此处为唯一放置）
    0);                 // 副本号

  // 在世界内继续构建细胞球（位于原点）
  DefineCell(world);

  // 打印世界体积的几何与材料信息
  PrintParameters(world);

  return world;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume*
DetectorConstruction::DefineCell(G4VPhysicalVolume* mother)
{
  /// 构建细胞体积（细胞膜边界球），并在其内部构建细胞核。
  /// @param mother 母体积（世界体积）
  /// @return 细胞物理体积指针

  // 细胞(膜边界)球：材料为水（胞质近似），中心位于原点
  // —— 细胞几何体（球体）——
  auto* solidCell = new G4Orb("Cell", fCellRadius);

  // —— 细胞逻辑体积 ——
  auto* logicalCell =
    new G4LogicalVolume(solidCell, fMat, "Cell", nullptr, nullptr, nullptr);

  // —— 细胞物理体积（置于母体积原点，无旋转）——
  G4VPhysicalVolume* cell = new G4PVPlacement(
    nullptr,                       // 母体积旋转矩阵
    G4ThreeVector(),               // 相对母体积的平移（原点）
    logicalCell,                   // 关联的逻辑体积
    "Cell",                        // 体积名称
    mother->GetLogicalVolume(),    // 母逻辑体积（世界）
    false,                         // 是否复制体积
    0);                            // 副本号

  // 在细胞内继续构建细胞核球（位于原点）
  DefineNucleus(cell);

  // 打印细胞体积的几何与材料信息
  PrintParameters(cell);

  return cell;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume*
DetectorConstruction::DefineNucleus(G4VPhysicalVolume* mother)
{
  /// 构建细胞核体积（球），敏感探测器将挂载于此体积。
  /// @param mother 母体积（细胞体积）
  /// @return 细胞核物理体积指针

  // 细胞核球：材料为水（核成分暂以水近似），中心位于原点
  // —— 细胞核几何体（球体）——
  auto* solidNucleus = new G4Orb("Nucleus", fNucleusRadius);

  // —— 细胞核逻辑体积 ——
  auto* logicalNucleus =
    new G4LogicalVolume(solidNucleus, fMat, "Nucleus", nullptr, nullptr, nullptr);

  // —— 细胞核物理体积（置于母体积原点，无旋转）——
  G4VPhysicalVolume* nucleus = new G4PVPlacement(
    nullptr,                       // 母体积旋转矩阵
    G4ThreeVector(),               // 相对母体积的平移（原点）
    logicalNucleus,                // 关联的逻辑体积
    "Nucleus",                     // 体积名称
    mother->GetLogicalVolume(),    // 母逻辑体积（细胞）
    false,                         // 是否复制体积
    0);                            // 副本号

  // 打印细胞核体积的几何与材料信息
  PrintParameters(nucleus);

  return nucleus;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::ConstructSDandField()
{
  /// 创建并注册敏感探测器(SD)，并将其挂载到细胞核逻辑体积上。
  /// 域级采样与核级打分均在细胞核内进行。

  // 敏感探测器名称与命中集合名称
  G4String sdName = "Nucleus";              // 探测器名称（对应细胞核体积）
  G4String hitCollName = "TrackerHitColl";  // 命中集合名称

  // 创建敏感探测器实例
  auto* trackerSD = new TrackerSD(sdName, hitCollName);

  // 向 SD 管理器注册该探测器
  G4SDManager* sdManager = G4SDManager::GetSDMpointer();
  sdManager->AddNewDetector(trackerSD);

  // 将探测器挂载到细胞核逻辑体积
  // 第三个参数 true 表示同时挂载到该体积的所有子体积
  SetSensitiveDetector("Nucleus", trackerSD, true);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetMaterial(const G4String& name)
{
  /// 按名称设置介质材料（通常由 DetectorMessenger 的命令调用）。
  /// 若指定材料不存在，则回退到 G4_WATER 并给出警告。
  /// @param name NIST 材料名称（如 "G4_WATER"）

  // 获取 NIST 材料管理器单例
  G4NistManager* nist = G4NistManager::Instance();

  // 按名称查找材料
  G4Material* material = nist->FindOrBuildMaterial(name);

  // 未找到指定材料：回退到 G4_WATER 并发出警告
  if (!material) {
    G4ExceptionDescription ed;
    ed << "未找到材料 '" << name << "'，回退使用 G4_WATER。";
    G4Exception("DetectorConstruction::SetMaterial", "MTK001", JustWarning, ed);
    material = nist->FindOrBuildMaterial("G4_WATER");
  }

  // 仅当材料确实改变时才更新（避免无意义的重复赋值）
  if (fMat != material) {
    fMat = material;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::PrintParameters(G4VPhysicalVolume* physVol) const
{
  /// 打印给定物理体积的几何与材料信息（名称、位置、形状、尺寸、材料）。
  /// @param physVol 待打印信息的物理体积指针

  // 空指针保护：避免解引用空指针导致崩溃
  if (!physVol) {
    G4cerr << "错误：物理体积指针为空！" << G4endl;
    return;
  }

  // 提取该体积的位置、几何实体与材料信息
  G4ThreeVector position = physVol->GetObjectTranslation();       // 体积中心位置
  G4LogicalVolume* logicalVolume = physVol->GetLogicalVolume();   // 关联逻辑体积
  G4VSolid* solid = logicalVolume->GetSolid();                    // 几何实体
  G4Material* material = logicalVolume->GetMaterial();            // 材料

  G4cout << "\n================ 体积参数 ================" << G4endl;
  G4cout << "物理体积名称: " << physVol->GetName() << G4endl;
  G4cout << "位置 (mm): " << position / mm << G4endl;
  G4cout << "材料: " << material->GetName() << G4endl;

  // 根据几何体类型，打印对应的尺寸信息
  if (auto* box = dynamic_cast<G4Box*>(solid)) {
    // 长方体：打印三个方向的完整边长
    G4double fullX = 2 * box->GetXHalfLength() / mm;  // x 方向完整边长
    G4double fullY = 2 * box->GetYHalfLength() / mm;  // y 方向完整边长
    G4double fullZ = 2 * box->GetZHalfLength() / mm;  // z 方向完整边长
    G4cout << "形状: 长方体" << G4endl;
    G4cout << "尺寸 (完整边长, mm): " << fullX << " x " << fullY << " x " << fullZ
           << G4endl;
  }
  else if (auto* orb = dynamic_cast<G4Orb*>(solid)) {
    // 球体：打印半径
    G4cout << "形状: 球体 (Orb)" << G4endl;
    G4cout << "半径 (um): " << orb->GetRadius() / um << G4endl;
  }
  else {
    // 未知或暂未支持的几何体类型
    G4cout << "形状: 未知或暂未支持的类型。" << G4endl;
  }

  G4cout << "===================================================\n" << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
