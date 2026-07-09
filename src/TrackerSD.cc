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
/// \file TrackerSD.cc
/// \brief TrackerSD 类的实现：细胞核敏感探测器
///
/// 挂载于细胞核逻辑体积，收集核内能量沉积命中(Hit)，并进行微剂量学单事件打分：
/// 域(site)随机抽样、线能 y、比能 z、核级比能 z_n 等，
/// 结果写入分析管理器的直方图与 ntuple，供后处理计算 S 值、LQ 参数等。

#include "TrackerSD.hh"

#include "G4AnalysisManager.hh"
#include "G4Event.hh"
#include "G4Exception.hh"
#include "G4ParticleDefinition.hh"
#include "G4PrimaryParticle.hh"
#include "G4PrimaryVertex.hh"
#include "G4ProcessType.hh"
#include "G4RunManager.hh"
#include "G4SDManager.hh"
#include "G4Step.hh"
#include "G4StepPoint.hh"
#include "G4SystemOfUnits.hh"
#include "G4VProcess.hh"
#include "Randomize.hh"
#include "globals.hh"

#include "DetectorConstruction.hh"
#include "EventAction.hh"
#include "PrimaryGeneratorAction.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

TrackerSD::TrackerSD(const G4String& sdName, const G4String& hitsCollectionName,
                     const int /*depthIndex*/)
  : G4VSensitiveDetector(sdName)
{
  /// 构造函数：设置敏感探测器名称并登记命中集合名称。
  /// @param sdName 敏感探测器名称（对应细胞核体积名称）
  /// @param hitsCollectionName 命中集合名称（注册到 collectionName 容器）
  /// @param depthIndex 层次深度参数（当前未使用）

  // 将命中集合名称登记到基类的 collectionName 容器
  collectionName.insert(hitsCollectionName);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

TrackerSD::~TrackerSD() = default;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void TrackerSD::Initialize(G4HCofThisEvent* hce)
{
  /// 事件开始时创建并注册本事件的命中集合。
  /// @param hce 本事件的命中集合容器（HCofThisEvent）

  // —— 创建命中集合 ——
  fHitsCollection =
    new TrackerHitColl(SensitiveDetectorName, collectionName[0]);

  // —— 向全局命中容器注册本集合 ——
  G4int collID =
    G4SDManager::GetSDMpointer()->GetCollectionID(collectionName[0]);

  hce->AddHitsCollection(collID, fHitsCollection);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool TrackerSD::ProcessHits(G4Step* aStep, G4TouchableHistory*)
{
  /// 每步步长处理：若有能量沉积则记录一个命中。
  /// 同时对初级粒子(parentID==0)记录核边界入射/出射动能(诊断量)。
  /// @param aStep 当前步对象
  /// @return 该步有非零能量沉积时返回 true，否则 false

  // —— 获取径迹与步端点 ——
  G4Track* aTrack = aStep->GetTrack();
  G4StepPoint* preStepPoint = aStep->GetPreStepPoint();    // 步前点
  G4StepPoint* postStepPoint = aStep->GetPostStepPoint();  // 步后点

  // 获取母粒子 ID，用于判定是否为初级粒子
  G4int parentID = aTrack->GetParentID();

  // —— 仅对初级粒子记录核边界入射/出射动能 ——
  // 注: 用 fGeomBoundary 状态位判定。对 DNA 物理下的 α 离子边界步该状态位可能不触发,
  //     故 α 的 T_in/T_out 可能为空(诊断量, 不影响 z/y/ε 打分)。
  //     若任务6.2 需精确能量平衡, 改用 G4UserSteppingAction 按体积过渡判定。
  if (parentID == 0) {
    if (preStepPoint->GetStepStatus() == fGeomBoundary) {
      G4AnalysisManager::Instance()->FillH1(12, preStepPoint->GetKineticEnergy());
    }
    if (postStepPoint->GetStepStatus() == fGeomBoundary) {
      G4AnalysisManager::Instance()->FillH1(13, postStepPoint->GetKineticEnergy());
    }
  }

  // —— 获取本步能量沉积；为零则跳过 ——
  G4double edep = aStep->GetTotalEnergyDeposit();
  if (edep == 0.) return false;

  // —— 创建新命中并填充能量与位置 ——
  auto newHit = new TrackerHit();
  newHit->SetEdep(edep);

  G4ThreeVector hitPosition = postStepPoint->GetPosition();  // 命中位置(步后点)
  newHit->SetPosition(hitPosition);

  // 将命中加入集合
  fHitsCollection->insert(newHit);

  return true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void TrackerSD::EndOfEvent(G4HCofThisEvent*)
{
  /// 事件结束时的微剂量学打分。
  /// 完成域(site)随机抽样、线能 y / 比能 z / 核级比能 z_n 计算，
  /// 填充直方图与逐事件 ntuple（hit/miss 事件均记录一行）。

  // —— 统计本事件命中数 ——
  const auto nofHits = fHitsCollection->entries();

  if (verboseLevel > 0) {
    G4cout << "  本事件命中集合共有 " << nofHits << " 个命中。" << G4endl;
  }

  // 注：miss 事件(nofHits==0, 如膜/胞外源 α 未命中核)不再早返回——
  // ntuple 用 hitFlag=0 记录(任务4.2)；直方图 f_{n,1} 由各 if(nucleusEdep>0)/if(nHint>0) 自然过滤。

  // —— 获取探测器几何参数（核半径、域半径、材料密度）——
  const DetectorConstruction* detConstruction =
    static_cast<const DetectorConstruction*>(
      G4RunManager::GetRunManager()->GetUserDetectorConstruction());

  G4double nucleusRadius = 0.;      // 细胞核半径
  G4double siteRadius = 0.;         // 域(site)半径
  G4double detectorDensity = 0.;    // 探测器材料密度
  G4ThreeVector hitPos;             // 选中的命中位置(用于 site 抽样)
  G4ThreeVector randCenterPos;      // site 球随机中心位置

  if (detConstruction) {
    siteRadius = detConstruction->GetSiteRadius();
    nucleusRadius = detConstruction->GetNucleusRadius();
    auto mat = G4Material::GetMaterial(detConstruction->GetMaterial());
    if (mat)
      detectorDensity = mat->GetDensity();
    else {
      G4ExceptionDescription msg;
      msg << "未找到材料，发生意外情况。" << G4endl;
      G4Exception("TrackerSD::EndOfEvent()", "TrackerSD002", JustWarning, msg);
    }
  }
  else {
    G4ExceptionDescription msg;
    msg << "未找到探测器构造，发生意外情况。" << G4endl;
    G4Exception("TrackerSD::EndOfEvent()", "TrackerSD001", JustWarning, msg);
  }

  // —— 事件级计数器与能量累加器 ——
  G4int nHsel = 0;            // 选择区(核内)命中数
  G4int nHsite = 0;           // 域(site)内命中数
  G4int nHint = 0;            // 同时在域内与选择区内的命中数
  G4double evtEdep = 0.;      // 域内能量沉积(单事件)
  G4double nucleusEdep = 0.;  // 核内总能量沉积(单事件, 任务4.1 用于 z_n)

  // —— 若有命中，随机选取一个命中作为 site 抽样锚点 ——
  if (nofHits > 0) {
    const G4int maxTries = 1000;  // 最大尝试次数
    G4int tries = 0;
    G4bool found = false;

    while (tries < maxTries && !found) {
      std::size_t randHit = static_cast<std::size_t>(G4UniformRand() * nofHits);
      auto hit = (*fHitsCollection)[randHit];
      G4ThreeVector hitPosition = hit->GetPosition();

      // 判定该 hit 是否落在细胞核内(核球为 hit 选择区)
      //   且要求 site 球能完全置于核内: 即此 hit 周围 r_d 球内的 site 中心
      //   满足 |randCenterPos| + r_d ≤ R_n. 由于 randCenterPos = hitPos + 偏移(|偏移|≤r_d),
      //   充要条件是 |hitPos| ≤ R_n - r_d. 这样保证 site 完全在核内,
      //   避免 site 部分超出核导致 nHint 偏小、nHsel/nHint 加权比偏高.
      if (hitPosition.mag() <= nucleusRadius - siteRadius)
      {
        hitPos = hitPosition;  // 记录有效的命中位置

        found = true;
      }
      ++tries;
    }

    if (found) {
      G4double siteRadius2 = siteRadius * siteRadius;  // 域半径平方

      // 在命中点周围随机放置 site 球心（球内均匀采样偏移量）
      G4double xRand, yRand, zRand, randRad2;
      do {
        xRand = (2 * G4UniformRand() - 1) * siteRadius;
        yRand = (2 * G4UniformRand() - 1) * siteRadius;
        zRand = (2 * G4UniformRand() - 1) * siteRadius;
        randRad2 = xRand * xRand + yRand * yRand + zRand * zRand;
      } while (randRad2 > siteRadius2);

      randCenterPos = G4ThreeVector(xRand + hitPos.x(), yRand + hitPos.y(),
                                    zRand + hitPos.z());
    }

    else {
      G4cout
        << "本事件尝试 " << tries << " 次后仍未在核内找到命中。" << G4endl
        << "请检查源位置 / 核半径是否合理。" << G4endl;
      // 不再 return: miss 事件仍需填 ntuple(hitFlag=0), 保持任务 4.2 hit/miss 都记录的设计
      // 此时 nHsel=nHsite=nHint=0, evtEdep=0; z_n/z_d 直方图被 if(nucleusEdep>0)/if(nHint>0) 自然过滤
    }
  }

  // —— 遍历所有命中，累加能量并统计各区域命中数 ——
  for (std::size_t jj = 0; jj < nofHits; jj++) {
    auto currentHit = (*fHitsCollection)[jj];
    G4ThreeVector currentHitPosition = currentHit->GetPosition();

    // 核总能量累加(所有 hit 均在核内, SD 在核上)
    nucleusEdep += currentHit->GetEdep();

    // 判定该 hit 是否落在细胞核内
    G4bool inNucleus = (currentHitPosition.mag() < nucleusRadius);

    // 判定该 hit 是否落在域(site)内
    G4double dist = (currentHitPosition - randCenterPos).mag();
    if (dist < siteRadius) {
      nHsite++;
      evtEdep += currentHit->GetEdep();
    }

    // 核内命中计数
    if (inNucleus) nHsel++;

    // 同时在域内且在核内的命中计数
    if (dist < siteRadius && inNucleus) nHint++;
  }

  // —— 获取分析管理器 ——
  auto analysisManager = G4AnalysisManager::Instance();

  // ===== P0 修复 #1: 补加出射边界步核内 edep =====
  //   SD hit-sum 已包含 pre/post 均在核内的 step + 入射 step (pre 在核外 post 在核内,
  //   edep 归到 post → SD 收到). 但出射 step (pre 在核内 post 在核外) 的 edep 归到
  //   post volume (核外), SD 不收到. SteppingAction 已通过 EventAction 累加了
  //   fNucleusEdepBoundary, 这里加到 nucleusEdep 得到完整核内沉积.
  {
    const auto* evtAct = dynamic_cast<const EventAction*>(
      G4RunManager::GetRunManager()->GetUserEventAction());
    if (evtAct) nucleusEdep += evtAct->GetNucleusEdepBoundary();
  }

  // ===== 核级比能 z_n (任务4.1)：核内总沉积 / 核质量，无加权(核是整个位点) =====
  G4double massNucleus =
    (4. / 3.) * CLHEP::pi * nucleusRadius * nucleusRadius * nucleusRadius
    * detectorDensity;        // 核质量
  G4double z_n = nucleusEdep / massNucleus;  // 核级比能
  if (nucleusEdep > 0.) {
    G4double zn = z_n / gray;
    G4int idF = analysisManager->GetH1Id("fzn");
    G4int idZF = analysisManager->GetH1Id("znfzn");
    G4int idZ2F = analysisManager->GetH1Id("z2nfzn");
    if (idF >= 0) analysisManager->FillH1(idF, zn, 1.);            // → z̄_{n,F}
    if (idZF >= 0) analysisManager->FillH1(idZF, zn, zn);          // → z̄_{n,D}
    if (idZ2F >= 0) analysisManager->FillH1(idZ2F, zn, zn * zn);   // 高阶矩
  }

  // —— 计算线能 y = ε / l̄，球体平均弦长 l̄ = 4r/3 ——
  G4double y = (evtEdep) / ((4. / 3.) * siteRadius);

  // —— 计算域级比能 z = ε / m_site ——
  G4double mass = ((4. / 3.) * CLHEP::pi * siteRadius * siteRadius
                   * siteRadius * detectorDensity);  // 域质量
  G4double z = (evtEdep / mass);                      // 域级比能

  // —— 填充微剂量学直方图（仅当域与核交集有命中时）——
  if (nHint > 0) {
    // 加权 = 选择区命中数 / 交集命中数
    G4double wght = G4double(nHsel) / G4double(nHint);

    // 直方图 0: 单事件授与能 (keV)
    analysisManager->FillH1(0, evtEdep / keV, wght);

    // 直方图 1: 加权单事件授与能 (keV)
    analysisManager->FillH1(1, evtEdep / keV, (evtEdep * wght / keV));

    // 直方图 2: 加权能量平方 (keV^2)
    analysisManager->FillH1(2, evtEdep / keV,
                            ((evtEdep * evtEdep * wght) / (keV * keV)));

    // 直方图 3: 线能 y (keV/um)
    analysisManager->FillH1(3, y / (keV / um), wght);

    // 直方图 4: 剂量加权线能 y (keV/um)
    analysisManager->FillH1(4, y / (keV / um), (y * wght / (keV / um)));

    // 直方图 5: 加权线能平方 ((keV/um)^2)
    analysisManager->FillH1(5, y / (keV / um),
                            (y * y * wght / ((keV / um) * (keV / um))));

    // 直方图 6: 比能 z (Gy)
    analysisManager->FillH1(6, z / gray, wght);

    // 直方图 7: 加权比能 z (Gy)
    analysisManager->FillH1(7, z / gray, (z * wght / gray));

    // 直方图 8: 加权比能平方 (Gy^2)
    analysisManager->FillH1(8, z / gray, (z * z * wght / (gray * gray)));

    // 直方图 9: 选择区命中数
    analysisManager->FillH1(9, nHsel, 1.);

    // 直方图 10: 域内命中数
    analysisManager->FillH1(10, nHsite, 1.);

    // 直方图 11: 域与选择区交集命中数
    analysisManager->FillH1(11, nHint, 1.);

    // 直方图 14: 域内命中数 vs 单事件授与能 (2D)
    analysisManager->FillH2(0, evtEdep / keV, nHsite, 1.);
  }
  // (else: miss 事件 nofHits==0 —— 对膜/胞外源属正常, 不再警告; ntuple 用 hitFlag=0 记录)

  // ===== 任务4.2：逐事件配对 ntuple (hit & miss 都填一行) =====
  // 取初级 α 动能(从 primary vertex)与区室编号(从 PrimaryGeneratorAction)
  G4double alphaE = 0.;
  const G4Event* evt = G4RunManager::GetRunManager()->GetCurrentEvent();
  if (evt && evt->GetNumberOfPrimaryVertex() > 0) {
    const G4PrimaryVertex* pv = evt->GetPrimaryVertex(0);
    if (pv) {
      const G4PrimaryParticle* pp = pv->GetPrimary(0);
      if (pp) alphaE = pp->GetKineticEnergy();
    }
  }

  // 区室编号映射
  G4int compId = 4;  // 未知区室
  const auto* pga = dynamic_cast<const PrimaryGeneratorAction*>(
    G4RunManager::GetRunManager()->GetUserPrimaryGeneratorAction());
  if (pga) {
    const G4String c = pga->GetCompartment();
    if (c == "Nucleus") compId = 0;
    else if (c == "Cytoplasm") compId = 1;
    else if (c == "Membrane") compId = 2;
    else if (c == "Extracellular") compId = 3;
  }

  // 加权（miss 事件取 0）
  G4double wght = (nHint > 0) ? G4double(nHsel) / G4double(nHint) : 0.;

  // —— 填充 ntuple 各列 ——
  analysisManager->FillNtupleIColumn(0, evt ? evt->GetEventID() : 0);
  analysisManager->FillNtupleDColumn(1, alphaE / MeV);
  analysisManager->FillNtupleDColumn(2, evtEdep / keV);
  analysisManager->FillNtupleDColumn(3, z / gray);          // z_d
  analysisManager->FillNtupleDColumn(4, nucleusEdep / keV);
  analysisManager->FillNtupleDColumn(5, z_n / gray);        // z_n
  analysisManager->FillNtupleDColumn(6, wght);
  analysisManager->FillNtupleIColumn(7, nHsel);
  analysisManager->FillNtupleIColumn(8, nHsite);
  analysisManager->FillNtupleIColumn(9, nHint);
  analysisManager->FillNtupleIColumn(10, nofHits > 0 ? 1 : 0);  // hitFlag
  analysisManager->FillNtupleIColumn(11, compId);               // compartment
  // 任务6.2: 全局能量沉积(从 EventAction 读, SteppingAction 每步累加)
  const auto* evtAct = dynamic_cast<const EventAction*>(
    G4RunManager::GetRunManager()->GetUserEventAction());
  G4double etot = evtAct ? evtAct->GetTotalEdep() : 0.;
  analysisManager->FillNtupleDColumn(12, etot / keV);

  // —— 提交本事件 ntuple 行 ——
  analysisManager->AddNtupleRow();
}
