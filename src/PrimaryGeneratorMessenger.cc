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
/// \file PrimaryGeneratorMessenger.cc
/// \brief PrimaryGeneratorMessenger 类的实现：注册初级粒子源相关的 UI 命令
///
/// 该 Messenger 向 UI 暴露 /beam/ 与 /source/ 目录下的命令，用于在运行时
/// 设置粒子源位置、源类型与分布区间，并将命令转发给 PrimaryGeneratorAction。

#include "PrimaryGeneratorMessenger.hh"

#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIdirectory.hh"

#include "PrimaryGeneratorAction.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorMessenger::PrimaryGeneratorMessenger(
  PrimaryGeneratorAction* primaryGen)
  : G4UImessenger(), fPrimaryAction(primaryGen)
{
  /// 构造函数：注册 /beam/ 与 /source/ 目录及其下的 UI 命令。
  /// 包括粒子源 Z 坐标、源类型(proton/ac225/alpha)以及 Ac-225 的分布区间。
  /// @param primaryGen 关联的 PrimaryGeneratorAction 对象指针，命令将转发给该对象

  // —— 束流(beam)命令目录 ——
  fGunDir = std::make_unique<G4UIdirectory>("/beam/");
  fGunDir->SetGuidance("初级粒子产生器控制命令");

  // —— 设置粒子源 Z 坐标命令 ——
  fZ0Cmd =
    std::make_unique<G4UIcmdWithADoubleAndUnit>("/beam/position/Z0", this);
  fZ0Cmd->SetGuidance("设置粒子源的 Z 坐标");
  fZ0Cmd->SetParameterName("posZ", false);
  fZ0Cmd->SetDefaultValue(0.0);
  fZ0Cmd->SetDefaultUnit("cm");
  fZ0Cmd->SetUnitCategory("Length");
  fZ0Cmd->AvailableForStates(G4State_PreInit, G4State_Idle);

  // —— 源(source)命令目录 ——
  fSourceDir = std::make_unique<G4UIdirectory>("/source/");
  fSourceDir->SetGuidance("粒子源类型控制命令");

  // —— 设置粒子源类型命令 ——
  fSourceTypeCmd = std::make_unique<G4UIcmdWithAString>("/source/type", this);
  fSourceTypeCmd->SetGuidance(
    "源类型：proton（基线）| ac225（Ac-225 alpha，分布区间由 /source/compartment 设定）"
    " | alpha（单能 alpha，能量由 /gun/energy 设定）"
    " | ac225_decay（路线2：静止 Ac-225 完整衰变链 4α+β+γ+反冲，任务7.1）"
    " | ac225_single_decay（路线2：单次 Ac-225 → Fr-221 + α(5.83 MeV)，α+反冲核同事件双顶点，任务X）"
    " | lu177_decay（路线2：静止 Lu-177 β⁻ + 208/113 keV γ + 反冲，任务8）");
  fSourceTypeCmd->SetParameterName("type", false);
  fSourceTypeCmd->SetCandidates("proton ac225 alpha ac225_decay ac225_single_decay lu177_decay");
  fSourceTypeCmd->SetDefaultValue("ac225");
  fSourceTypeCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

  // —— 设置 Ac-225 分布区间命令 ——
  fCompartmentCmd = std::make_unique<G4UIcmdWithAString>("/source/compartment", this);
  fCompartmentCmd->SetGuidance(
    "源分布区间(ac225): Nucleus | Cytoplasm | Membrane | Extracellular"
    " | WholeCell (整个细胞均匀 = 核∪质∪膜, 即 Rc 球内)"
    " | CellExceptNucleus (除核外均匀 = 质∪膜, 不含核)");
  fCompartmentCmd->SetParameterName("comp", false);
  fCompartmentCmd->SetCandidates("Nucleus Cytoplasm Membrane Extracellular WholeCell CellExceptNucleus");
  fCompartmentCmd->SetDefaultValue("Membrane");
  fCompartmentCmd->AvailableForStates(G4State_PreInit, G4State_Idle);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorMessenger::~PrimaryGeneratorMessenger() = default;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PrimaryGeneratorMessenger::SetNewValue(G4UIcommand* command,
                                            G4String newValue)
{
  /// 命令分发：根据 command 指针将 UI 传入的新值转发给 PrimaryGeneratorAction。
  /// 注意：此处使用独立 if（非 else if），保持与原逻辑一致。
  /// @param command 指向触发的 UI 命令对象
  /// @param newValue UI 传入的命令参数字符串

  // 粒子源 Z 坐标
  if (command == fZ0Cmd.get()) {
    fPrimaryAction->SetPositionZ(fZ0Cmd->GetNewDoubleValue(newValue));
  }
  // 粒子源类型
  if (command == fSourceTypeCmd.get()) {
    fPrimaryAction->SetSourceType(newValue);
  }
  // Ac-225 分布区间
  if (command == fCompartmentCmd.get()) {
    fPrimaryAction->SetCompartment(newValue);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
