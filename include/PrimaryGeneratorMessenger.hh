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
/// \file PrimaryGeneratorMessenger.hh
/// \brief PrimaryGeneratorMessenger 类的定义：粒子源的 UI 命令交互

#ifndef PrimaryGeneratorMessenger_h
#define PrimaryGeneratorMessenger_h 1

#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UImessenger.hh"
#include "globals.hh"


class PrimaryGeneratorAction;
class G4UIdirectory;
class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWithAString;

#include <memory>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

/// @brief 粒子源交互命令类
///
/// 定义 /beam/ 与 /source/ 两个命令目录：前者控制质子枪源点 Z 坐标，
/// 后者控制 Ac-225 α 源的类型与区室。命令转发给 PrimaryGeneratorAction。
class PrimaryGeneratorMessenger : public G4UImessenger
{
  public:
    PrimaryGeneratorMessenger(PrimaryGeneratorAction*);  // 注册 /beam/ 与 /source/ 命令
    ~PrimaryGeneratorMessenger() override;               // 析构（命令对象由智能指针释放）
    /// 响应 UI 命令赋值，转发给 PrimaryGeneratorAction 的对应方法。
    /// @param command 触发的 UI 命令指针
    /// @param newValue 用户输入的新值字符串
    void SetNewValue(G4UIcommand*, G4String) override;

  private:
    // ===== 命令目标与命令对象 =====
    PrimaryGeneratorAction* fPrimaryAction = nullptr;                  // 关联的粒子源动作对象
    std::unique_ptr<G4UIdirectory> fGunDir;                            // /beam/ 命令目录
    std::unique_ptr<G4UIdirectory> fSourceDir;                         // /source/ 命令目录
    std::unique_ptr<G4UIcmdWithADoubleAndUnit> fZ0Cmd;                 // /beam/position/Z0 命令
    std::unique_ptr<G4UIcmdWithAString> fSourceTypeCmd;                // /source/type 命令
    std::unique_ptr<G4UIcmdWithAString> fCompartmentCmd;               // /source/compartment 命令
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
