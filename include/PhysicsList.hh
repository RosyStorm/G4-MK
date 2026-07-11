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
/// \file PhysicsList.hh
/// \brief PhysicsList 类的定义：物理过程列表

#ifndef PhysicsList_h
#define PhysicsList_h 1

#include "G4SystemOfUnits.hh"
#include "G4VModularPhysicsList.hh"
#include "globals.hh"

#include <memory>

class PhysicsListMessenger;

//.....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

/// @brief 物理过程列表类
///
/// 基于 G4VModularPhysicsList，默认采用 G4EmDNAPhysics_option2（Geant4-DNA
/// option2）作为电磁物理构造器，并构造全套粒子（含 DNA 离子：hydrogen、
/// deuteron、triton、helium、alpha、alpha+、alpha++、carbon）。运行时可通过
/// PhysicsListMessenger 的 /physics/addPhysics 命令切换物理构造器。
class PhysicsList : public G4VModularPhysicsList
{
  public:
    // ===== 构造与析构 =====
    PhysicsList();   // 创建 Messenger、设置默认 DNA 物理(option2)
    ~PhysicsList() override;  // 析构（Messenger 与物理构造器由智能指针释放）

    // ===== Geant4 强制重载接口 =====
    void ConstructParticle() override;  // 构造全部粒子定义
    void ConstructProcess() override;   // 构造粒子输运与电磁物理过程

    // ===== 物理构造器切换 =====
    /// 按名称切换物理构造器（dna_opt1..8、liv、penelope、em_standard_opt4 等）。
    /// @param name 物理构造器名称字符串
    void AddPhysicsList(const G4String& name);

  private:
    // ===== 物理构造器与交互命令 =====
    G4String fName = "";                                       // 当前激活的物理构造器名称
    std::unique_ptr<G4VPhysicsConstructor> fPhysicsList;       // 当前电磁物理构造器
    std::unique_ptr<G4VPhysicsConstructor> fDecayPhysics;      // 放射性衰变物理(路线2, 任务7.1)
    std::unique_ptr<PhysicsListMessenger> fMessenger;          // UI 命令交互对象
};

//.....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

#endif
