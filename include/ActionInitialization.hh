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
/// \file ActionInitialization.hh
/// \brief ActionInitialization 类的定义：用户动作初始化

#ifndef ActionInitialization_h
#define ActionInitialization_h 1

#include "G4VUserActionInitialization.hh"

class DetectorConstruction;
class PhysicsList;

//.....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

/// @brief 用户动作初始化类
///
/// 负责为每个工作线程创建用户动作对象（PrimaryGeneratorAction、RunAction、
/// EventAction、SteppingAction）。在多线程运行模式下由 RunManager 针对每个
/// 工作线程调用 Build()，对主线程调用 BuildForMaster()。
class ActionInitialization : public G4VUserActionInitialization
{
  public:
    // ===== 构造与析构 =====
    ActionInitialization() = default;            // 默认构造
    ~ActionInitialization() override = default;  // 默认析构

    // ===== Geant4 强制重载接口 =====
    void BuildForMaster() const override;  // 主线程动作构建（仅注册 RunAction）
    void Build() const override;           // 工作线程动作构建（注册全部用户动作）
};

//.....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

#endif
