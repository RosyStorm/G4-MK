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
/// \file PhysicsListMessenger.hh
/// \brief PhysicsListMessenger 类的定义：物理列表的 UI 命令交互

#ifndef PhysicsListMessenger_h
#define PhysicsListMessenger_h 1

#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UImessenger.hh"
#include "globals.hh"


class PhysicsList;
class G4UIdirectory;

#include <memory>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

/// @brief 物理列表交互命令类
///
/// 定义 /physics/ 命令目录下的 UI 命令，使用户可在运行时（宏文件或交互式会话）
/// 切换电磁物理构造器。命令转发给 PhysicsList::AddPhysicsList。
class PhysicsListMessenger : public G4UImessenger
{
  public:
    PhysicsListMessenger(PhysicsList*);   // 注册 /physics/ 命令
    ~PhysicsListMessenger() override;     // 析构（命令对象由智能指针释放）
    /// 响应 UI 命令赋值，转发给 PhysicsList 的对应方法。
    /// @param command 触发的 UI 命令指针
    /// @param newVal 用户输入的新值字符串
    void SetNewValue(G4UIcommand*, G4String) override;

  private:
    // ===== 命令目标与命令对象 =====
    PhysicsList* fPL = nullptr;                              // 关联的物理列表对象
    std::unique_ptr<G4UIdirectory> fPhysDir;                 // /physics/ 命令目录
    std::unique_ptr<G4UIcmdWithAString> fPhysicsListCmd;     // /physics/addPhysics 命令
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
