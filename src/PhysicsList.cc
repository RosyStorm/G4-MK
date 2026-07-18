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
/// \file PhysicsList.cc
/// \brief PhysicsList 类的实现：物理过程列表
///
/// 基于 Geant4-DNA 模块化物理构造器，注册粒子与电磁/DNA 物理过程。
/// 可通过 UI 命令在多个 DNA 物理选项(dna_opt1~8)及 Livermore/Penelope/标准电磁之间切换。

#include "PhysicsList.hh"

#include "G4SystemOfUnits.hh"
#include "G4UnitsTable.hh"

#include "PhysicsListMessenger.hh"

// 方案C: 混合区域物理(空气标准EM + 细胞DNA)
#include "DetectorConstruction.hh"
#include "G4RunManager.hh"
#include "G4LossTableManager.hh"
#include "G4EmConfigurator.hh"
#include "G4RegionStore.hh"
#include "G4BetheBlochModel.hh"
#include "G4MollerBhabhaModel.hh"
#include "G4UrbanMscModel.hh"

// 粒子定义
#include "G4Gamma.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"
#include "G4Proton.hh"
#include "G4BaryonConstructor.hh"
#include "G4BosonConstructor.hh"
#include "G4DNAGenericIonsManager.hh"
#include "G4IonConstructor.hh"
#include "G4LeptonConstructor.hh"
#include "G4MesonConstructor.hh"
#include "G4ShortLivedConstructor.hh"

// 物理包（构建器包含在 Geant4 源码中）
// 电磁
#include "G4EmDNAPhysics.hh"
#include "G4EmDNAPhysics_option1.hh"
#include "G4EmDNAPhysics_option2.hh"
#include "G4EmDNAPhysics_option3.hh"
#include "G4EmDNAPhysics_option4.hh"
#include "G4EmDNAPhysics_option5.hh"
#include "G4EmDNAPhysics_option6.hh"
#include "G4EmDNAPhysics_option7.hh"
#include "G4EmDNAPhysics_option8.hh"
#include "G4EmLivermorePhysics.hh"
#include "G4EmPenelopePhysics.hh"
#include "G4EmStandardPhysics_option4.hh"
#include "G4RadioactiveDecayPhysics.hh"   // 放射性衰变(路线2, 任务7.1)



#include <memory>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PhysicsList::PhysicsList() : G4VModularPhysicsList()
{
  /// 构造函数：创建交互命令对象、设置默认物理包(dna_opt2)。

  // 创建物理列表的交互命令对象
  fMessenger = std::make_unique<PhysicsListMessenger>(this);

  // 设置冗余级别
  SetVerboseLevel(1);

  // 默认使用 DNA 物理 option2
  fPhysicsList = std::make_unique<G4EmDNAPhysics_option2>();

  // 放射性衰变物理(路线2, 任务7.1): 自动模拟 Ac-225 完整衰变链(α/β/γ/核反冲)。
  // 常驻注册——仅影响放射性核素(Ac-225 及子核), 稳定的 α/质子不受影响。
  fDecayPhysics = std::make_unique<G4RadioactiveDecayPhysics>();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PhysicsList::~PhysicsList() = default;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysicsList::ConstructParticle()
{
  /// 注册所有需要的粒子定义（Geant4 应用启动时调用，勿与 run 初始化混淆）。
  /// 依次构造玻色子、轻子、介子、重子、离子及 DNA 专用离子。

  // —— 玻色子（光子等）——
  G4BosonConstructor pBosonConstructor;
  pBosonConstructor.ConstructParticle();

  // —— 轻子（电子、μ 子等）——
  G4LeptonConstructor pLeptonConstructor;
  pLeptonConstructor.ConstructParticle();

  // —— 介子 ——
  G4MesonConstructor pMesonConstructor;
  pMesonConstructor.ConstructParticle();

  // —— 重子（质子、中子等）——
  G4BaryonConstructor pBaryonConstructor;
  pBaryonConstructor.ConstructParticle();

  // —— 离子 ——
  G4IonConstructor pIonConstructor;
  pIonConstructor.ConstructParticle();

  // —— 短寿命粒子 ——
  G4ShortLivedConstructor pShortLivedConstructor;
  pShortLivedConstructor.ConstructParticle();

  // —— DNA 专用离子管理器：注册氢/氘/氚/氦/α 及带电 α 等 ——
  G4DNAGenericIonsManager* genericIonsManager;
  genericIonsManager = G4DNAGenericIonsManager::Instance();
  genericIonsManager->GetIon("hydrogen");
  genericIonsManager->GetIon("deuteron");
  genericIonsManager->GetIon("triton");
  genericIonsManager->GetIon("helium");
  genericIonsManager->GetIon("alpha");
  genericIonsManager->GetIon("alpha+");
  genericIonsManager->GetIon("alpha++");
  genericIonsManager->GetIon("carbon");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysicsList::ConstructProcess()
{
  /// 构造所有物理过程：输运 + 电磁/DNA 过程 + 截断设置。
  /// 方案C(Am-241): 若 fWorldIsAir, 注册标准EM+DNA 双物理,
  ///   用 G4EmConfigurator 在 CellRegion 中关闭标准EM(避免双重计数)。

  // —— 粒子输运过程 ——
  AddTransportation();

  // —— 检测是否 Am-241 混合物理模式 ——
  const auto* det = static_cast<const DetectorConstruction*>(
    G4RunManager::GetRunManager()->GetUserDetectorConstruction());
  G4bool mixedPhysics = det && det->GetWorldIsAir();

  if (mixedPhysics) {
    // ===== 方案C: 空气用标准EM, 细胞用DNA =====
    G4cout << "[PhysicsList] 混合物理模式: 空气=标准EM(option4), 细胞=DNA(option2)" << G4endl;

    // (a) 标准 EM(空气中 α 能损由 Bethe-Bloch 处理)
    fEmPhysics = std::make_unique<G4EmStandardPhysics_option4>();
    fEmPhysics->ConstructProcess();

    // (b) DNA 物理(细胞水中逐电离跟踪)
    fPhysicsList->ConstructProcess();

    // (c) 在 CellRegion 中关闭标准 EM 模型(避免水中的双重计数)
    //   仅当 CellRegion 存在时执行(master 在 Construct() 中创建, worker 通过共享几何可见)
    G4Region* cellRegion = G4RegionStore::GetInstance()->GetRegion("CellRegion");
    if (cellRegion) {
      auto* emConf = G4LossTableManager::Instance()->EmConfigurator();
      const G4double thresh = 9.9 * MeV;

      // alpha: ionIoni(标准) 在 CellRegion 关闭
      {
        auto* mod = new G4BetheBlochModel();
        mod->SetActivationLowEnergyLimit(thresh);
        emConf->SetExtraEmModel("alpha", "ionIoni", mod, "CellRegion");
      }
      // e-: eIoni + msc 在 CellRegion 关闭
      {
        auto* mod = new G4MollerBhabhaModel();
        mod->SetActivationLowEnergyLimit(thresh);
        emConf->SetExtraEmModel("e-", "eIoni", mod, "CellRegion");
      }
      {
        auto* mod = new G4UrbanMscModel();
        mod->SetActivationLowEnergyLimit(thresh);
        emConf->SetExtraEmModel("e-", "msc", mod, "CellRegion");
      }
      // proton: hIoni 在 CellRegion 关闭
      {
        auto* mod = new G4BetheBlochModel();
        mod->SetActivationLowEnergyLimit(thresh);
        emConf->SetExtraEmModel("proton", "hIoni", mod, "CellRegion");
      }

      emConf->AddModels();
      G4cout << "[PhysicsList] CellRegion: 标准 EM 已关闭, DNA 生效" << G4endl;
    }
  }
  else {
    // ===== 其他源: 纯 DNA(不变) =====
    fPhysicsList->ConstructProcess();
  }

  // —— 放射性衰变(路线2, 任务7.1): α/β/γ/核反冲衰变链 ——
  fDecayPhysics->ConstructProcess();

  // —— 设置粒子产生阈截断(cuts) ——
  SetCuts();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysicsList::AddPhysicsList(const G4String& name)
{
  /// 按名称切换物理构造器。
  /// @param name 物理列表名称（如 "dna_opt2"、"liv"、"penelope"、"em_standard_opt4" 等）

  if (verboseLevel > -1) {
    G4cout << "PhysicsList::AddPhysicsList: <" << name << ">" << G4endl;
  }

  // 与当前物理列表相同则不切换
  if (name == fName) return;

  // —— DNA 物理 ——
  if (name == "dna") {
    fName = name;
    fPhysicsList = std::make_unique<G4EmDNAPhysics>();
    G4cout << fPhysicsList->GetPhysicsName()
           << " 物理包已激活。" << G4endl;
  }
  else if (name == "dna_opt1") {
    fName = name;
    fPhysicsList = std::make_unique<G4EmDNAPhysics_option1>();
    G4cout << fPhysicsList->GetPhysicsName()
           << " 物理包已激活。" << G4endl;
  }

  else if (name == "dna_opt2") {
    fName = name;
    fPhysicsList = std::make_unique<G4EmDNAPhysics_option2>();
    G4cout << fPhysicsList->GetPhysicsName()
           << " 物理包已激活。" << G4endl;
  }

  else if (name == "dna_opt3") {
    fName = name;
    fPhysicsList = std::make_unique<G4EmDNAPhysics_option3>();
    G4cout << fPhysicsList->GetPhysicsName()
           << " 物理包已激活。" << G4endl;
  }

  else if (name == "dna_opt4") {
    fName = name;
    fPhysicsList = std::make_unique<G4EmDNAPhysics_option4>();
    G4cout << fPhysicsList->GetPhysicsName()
           << " 物理包已激活。" << G4endl;
  }

  else if (name == "dna_opt5") {
    fName = name;
    fPhysicsList = std::make_unique<G4EmDNAPhysics_option5>();
    G4cout << fPhysicsList->GetPhysicsName()
           << " 物理包已激活。" << G4endl;
  }

  else if (name == "dna_opt6") {
    fName = name;
    fPhysicsList = std::make_unique<G4EmDNAPhysics_option6>();
    G4cout << fPhysicsList->GetPhysicsName()
           << " 物理包已激活。" << G4endl;
  }

  else if (name == "dna_opt7") {
    fName = name;
    fPhysicsList = std::make_unique<G4EmDNAPhysics_option7>();
    G4cout << fPhysicsList->GetPhysicsName()
           << " 物理包已激活。" << G4endl;
  }

  else if (name == "dna_opt8") {
    fName = name;
    fPhysicsList = std::make_unique<G4EmDNAPhysics_option8>();
    G4cout << fPhysicsList->GetPhysicsName()
           << " 物理包已激活。" << G4endl;
  }

  // —— Livermore 物理 ——
  else if (name == "liv") {
    fName = name;
    fPhysicsList = std::make_unique<G4EmLivermorePhysics>();
    G4cout << fPhysicsList->GetPhysicsName()
           << " 物理包已激活。" << G4endl;
  }

  // —— Penelope 物理 ——
  else if (name == "penelope") {
    fName = name;
    fPhysicsList = std::make_unique<G4EmPenelopePhysics>();
    G4cout << fPhysicsList->GetPhysicsName()
           << " 物理包已激活。" << G4endl;
  }

  // —— 标准电磁物理 option4 ——
  else if (name == "em_standard_opt4") {
    fName = name;
    fPhysicsList = std::make_unique<G4EmStandardPhysics_option4>();
    G4cout << fPhysicsList->GetPhysicsName()
           << " 物理包已激活。" << G4endl;
  }

  else {
    G4cout << "PhysicsList::AddPhysicsList: \"" << name << "\" 未定义!"
           << G4endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
