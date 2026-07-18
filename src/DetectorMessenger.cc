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
/// \file DetectorMessenger.cc
/// \brief DetectorMessenger 类的实现：注册探测器几何相关的 UI 命令
///
/// 该 Messenger 向 UI 暴露 /mygeom/ 目录下的命令，用于在运行时修改
/// 探测器几何参数（次级电子射程、各半径、材料）以及粒子杀死策略，
/// 并将命令转发给 DetectorConstruction。

#include "DetectorMessenger.hh"

#include "G4UIcmdWithABool.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithAString.hh"
#include "globals.hh"

#include "DetectorConstruction.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorMessenger::DetectorMessenger(DetectorConstruction* myDet)
  : fDetector(myDet)
{
  /// 构造函数：注册 /mygeom/ 目录及其下的全部 UI 命令。
  /// 包括次级电子最大射程、细胞半径、细胞核半径、域半径、材料名称，
  /// 以及细胞外杀死开关与核边界杀死开关。
  /// @param myDet 关联的 DetectorConstruction 对象指针，命令将转发给该对象

  // —— 探测器命令目录 ——
  fDetectorDir = std::make_unique<G4UIdirectory>("/mygeom/");
  fDetectorDir->SetGuidance("探测器控制命令。");

  // —— 设置次级电子最大射程命令 ——
  fMaxRangeCmd =
    std::make_unique<G4UIcmdWithADoubleAndUnit>("/mygeom/maxRange", this);
  fMaxRangeCmd->SetGuidance("次级电子的最大射程。");
  fMaxRangeCmd->SetGuidance("该参数用于设定");
  fMaxRangeCmd->SetGuidance("世界体积的尺寸。");
  fMaxRangeCmd->SetParameterName("range", false);
  fMaxRangeCmd->SetDefaultUnit("mm");
  fMaxRangeCmd->SetUnitCandidates("um mm cm");
  fMaxRangeCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

  // —— 设置细胞半径命令 ——
  fCellRadiusCmd =
    std::make_unique<G4UIcmdWithADoubleAndUnit>("/mygeom/cellRadius", this);
  fCellRadiusCmd->SetGuidance("细胞（细胞膜边界）的半径");
  fCellRadiusCmd->SetParameterName("radius", false);
  fCellRadiusCmd->SetDefaultUnit("um");
  fCellRadiusCmd->SetUnitCandidates("nm um mm cm");
  fCellRadiusCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

  // —— 设置介质材料名称命令 ——
  fMatNameCmd = std::make_unique<G4UIcmdWithAString>("/mygeom/material", this);
  fMatNameCmd->SetGuidance("介质材料的名称。");
  fMatNameCmd->SetParameterName("name", false);
  fMatNameCmd->AvailableForStates(G4State_PreInit);

  // —— 设置细胞核半径命令 ——
  fNucleusRadiusCmd =
    std::make_unique<G4UIcmdWithADoubleAndUnit>("/mygeom/nucleusRadius", this);
  fNucleusRadiusCmd->SetGuidance("细胞核（域打分区域）的半径");
  fNucleusRadiusCmd->SetParameterName("radius", false);
  fNucleusRadiusCmd->SetDefaultUnit("um");
  fNucleusRadiusCmd->SetUnitCandidates("nm um mm cm");
  fNucleusRadiusCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

  // —— 设置域(site)半径命令 ——
  fSiteRadiusCmd =
    std::make_unique<G4UIcmdWithADoubleAndUnit>("/mygeom/siteRadius", this);
  fSiteRadiusCmd->SetGuidance("域(site)的半径。");
  fSiteRadiusCmd->SetParameterName("radius", false);
  fSiteRadiusCmd->SetDefaultUnit("um");
  fSiteRadiusCmd->SetUnitCandidates("nm um mm cm");
  fSiteRadiusCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

  // —— 细胞外粒子杀死开关命令 ——
  fKillOutsideCellCmd =
    std::make_unique<G4UIcmdWithABool>("/mygeom/killOutsideCell", this);
  fKillOutsideCellCmd->SetGuidance(
    "加速开关：粒子一旦离开细胞(R_cell)且向外飞行即予以杀死。"
    " 生产运行默认 true；进行 alpha 射程验证(任务 3.1)时设为 false。");
  fKillOutsideCellCmd->SetParameterName("flag", false);
  fKillOutsideCellCmd->SetDefaultValue(true);
  fKillOutsideCellCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

  // —— 核边界杀死半径开关命令 ——
  fKillAtNucleusCmd =
    std::make_unique<G4UIcmdWithABool>("/mygeom/killAtNucleus", this);
  fKillAtNucleusCmd->SetGuidance(
    "杀死半径：true=按细胞核半径 R_n（默认，更快，不影响核内打分），"
    " false=按细胞半径 R_cell。");
  fKillAtNucleusCmd->SetParameterName("flag2", false);
  fKillAtNucleusCmd->SetDefaultValue(true);
  fKillAtNucleusCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

  // —— 世界材料空气开关 (Am-241 外照射 SSD 空气间隙) ——
  fWorldAirCmd =
    std::make_unique<G4UIcmdWithABool>("/mygeom/worldAir", this);
  fWorldAirCmd->SetGuidance(
    "世界材料: true=空气(Am-241 外照射, 细胞外为空气), false=水(默认, 其他源不变)。"
    " 细胞和细胞核始终为水。");
  fWorldAirCmd->SetParameterName("air", false);
  fWorldAirCmd->SetDefaultValue(false);
  fWorldAirCmd->AvailableForStates(G4State_PreInit);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorMessenger::~DetectorMessenger() = default;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorMessenger::SetNewValue(G4UIcommand* command, G4String newValue)
{
  /// 命令分发：根据 command 指针将 UI 传入的新值转发给 DetectorConstruction。
  /// @param command 指向触发的 UI 命令对象
  /// @param newValue UI 传入的命令参数字符串

  // 次级电子最大射程
  if (command == fMaxRangeCmd.get()) {
    fDetector->SetMaxRange(fMaxRangeCmd->GetNewDoubleValue(newValue));
  }
  // 介质材料名称
  else if (command == fMatNameCmd.get()) {
    fDetector->SetMaterial(newValue);
  }
  // 细胞半径
  else if (command == fCellRadiusCmd.get()) {
    fDetector->SetCellRadius(fCellRadiusCmd->GetNewDoubleValue(newValue));
  }
  // 细胞核半径
  else if (command == fNucleusRadiusCmd.get()) {
    fDetector->SetNucleusRadius(fNucleusRadiusCmd->GetNewDoubleValue(newValue));
  }
  // 域(site)半径
  else if (command == fSiteRadiusCmd.get()) {
    fDetector->SetSiteRadius(fSiteRadiusCmd->GetNewDoubleValue(newValue));
  }
  // 细胞外粒子杀死开关
  else if (command == fKillOutsideCellCmd.get()) {
    fDetector->SetKillOutsideCell(fKillOutsideCellCmd->GetNewBoolValue(newValue));
  }
  // 核边界杀死半径开关
  else if (command == fKillAtNucleusCmd.get()) {
    fDetector->SetKillAtNucleus(fKillAtNucleusCmd->GetNewBoolValue(newValue));
  }
  // 世界材料空气开关
  else if (command == fWorldAirCmd.get()) {
    fDetector->SetWorldIsAir(fWorldAirCmd->GetNewBoolValue(newValue));
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
