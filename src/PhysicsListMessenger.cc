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
/// \file PhysicsListMessenger.cc
/// \brief PhysicsListMessenger 类的实现：注册物理列表相关的 UI 命令
///
/// 该 Messenger 向 UI 暴露 /physics/ 目录下的命令，用于在运行时
/// 选择/切换物理模型，并将命令转发给 PhysicsList。

#include "PhysicsListMessenger.hh"

#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIdirectory.hh"

#include "PhysicsList.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PhysicsListMessenger::PhysicsListMessenger(PhysicsList* phys)
  : G4UImessenger(), fPL(phys)
{
  /// 构造函数：注册 /physics/ 目录及其下的 UI 命令。
  /// 当前注册物理列表选择命令(addPhysics)，用于按名称添加物理模型。
  /// @param phys 关联的 PhysicsList 对象指针，命令将转发给该对象

  // —— 物理命令目录 ——
  fPhysDir = std::make_unique<G4UIdirectory>("/physics/");
  fPhysDir->SetGuidance("用于设置物理模型与产生截断的命令");

  // —— 添加物理列表命令 ——
  fPhysicsListCmd =
    std::make_unique<G4UIcmdWithAString>("/physics/addPhysics", this);
  fPhysicsListCmd->SetGuidance("按名称添加物理列表。");
  fPhysicsListCmd->SetParameterName("PList", false);
  fPhysicsListCmd->AvailableForStates(G4State_PreInit);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PhysicsListMessenger::~PhysicsListMessenger() = default;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysicsListMessenger::SetNewValue(G4UIcommand* command, G4String newVal)
{
  /// 命令分发：根据 command 指针将 UI 传入的新值转发给 PhysicsList。
  /// @param command 指向触发的 UI 命令对象
  /// @param newVal UI 传入的命令参数字符串

  // 物理列表名称
  if (command == fPhysicsListCmd.get()) {
    fPL->AddPhysicsList(newVal);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
