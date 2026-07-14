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
/// \file microtrack.cc
/// \brief microtrack 示例的主程序

// 以下为 Geant4-DNA 协作组的示例来源与论文引用声明（须原样保留）：
// 本示例由 Geant4-DNA 协作组提供，公开发表基于 Geant4-DNA 的结果时须引用以下论文。
// This example is provided by the Geant4-DNA collaboration
// Any report or published results obtained using the Geant4-DNA software
// shall cite the following Geant4-DNA collaboration publications:
// Med. Phys. 45, (2018) e722-e739
// Phys. Med. 31 (2015) 861-874
// Med. Phys. 37 (2010) 4692-4708
// Int. J. Model. Simul. Sci. Comput. 1 (2010) 157–178
//
// Contact author: M.A. Cortes-Giraldo, Universidad de Sevilla, Spain
//                 (miancortes@us.es)
// Other authors: M. Galocha-Oliva, A. Baratto-Roldan
//
// Example ref. paper: A. Baratto-Roldan et al., Front Phys 9 (2021) 726787
//

#include "G4RunManagerFactory.hh"
#include "G4SteppingVerbose.hh"
#include "G4Types.hh"
#include "G4UIExecutive.hh"
#include "G4UImanager.hh"
#include "G4VisExecutive.hh"
#include "Randomize.hh"

#include "ActionInitialization.hh"
#include "DetectorConstruction.hh"
#include "PhysicsList.hh"
#include "RunAction.hh"

#include <fstream>
#include <regex>
#include <string>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

/// 从宏文件里抽出 /source/type xxx 与 /source/compartment yyy 设置。
/// master 线程不能用 GetUserPrimaryGeneratorAction()（MT 模式下 PGA 在 worker 上），
/// 必须从主程序拿 argv[1] 自己解析。兜底与 PrimaryGeneratorAction 默认一致。
/// @param macroPath 宏文件路径
static void ScanMacroForRunMeta(const std::string& macroPath)
{
  std::ifstream in(macroPath);
  if (!in.is_open()) {
    G4cout << "[RunMeta] 警告: 无法打开宏文件 " << macroPath
           << " — 使用默认源类型 ac225_decay / 区室 Membrane" << G4endl;
    return;
  }

  G4String srcType      = "ac225_decay";
  G4String compartment  = "Membrane";
  std::string line;
  // /source/type <xxx>          xxx ∈ {proton, ac225, alpha, ac225_decay, lu177_decay}
  // /source/compartment <yyy>   yyy ∈ {Nucleus, Cytoplasm, Membrane, Extracellular}
  std::regex reType(R"(^\s*/source/type\s+(\S+))");
  std::regex reComp(R"(^\s*/source/compartment\s+(\S+))");
  std::smatch m;

  while (std::getline(in, line)) {
    // 跳过注释
    if (line.empty() || line[0] == '#') continue;
    if (std::regex_search(line, m, reType) && m.size() >= 2) {
      srcType = m[1].str().c_str();
    }
    if (std::regex_search(line, m, reComp) && m.size() >= 2) {
      compartment = m[1].str().c_str();
    }
  }
  in.close();

  RunAction::SetRunMeta(srcType, compartment);
  G4cout << "[RunMeta] 从 " << macroPath
         << " 解析到: sourceType=" << srcType
         << ", compartment=" << compartment << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

/// 程序入口：初始化运行管理器、注册用户初始化类，并按交互/批处理模式启动模拟。
/// @param argc 命令行参数个数（1=交互模式，2=批处理，3=指定线程数）
/// @param argv 命令行参数数组（argv[1]=宏文件名，argv[2]=线程数）
/// @return 程序退出码
int main(int argc, char** argv)
{
  // —— 检测交互模式（无参数时）并创建 UI 会话 ——
  G4UIExecutive* ui = nullptr;
  if (argc == 1) {
    ui = new G4UIExecutive(argc, argv);
  }

  // —— 选择随机数引擎 ——
  G4Random::setTheEngine(new CLHEP::RanecuEngine);

  // —— 设置步进输出的精度并启用最佳单位显示 ——
  G4int precision = 4;
  G4SteppingVerbose::UseBestUnit(precision);

  // —— 创建运行管理器（默认类型，支持多线程）——
  auto runManager =
    G4RunManagerFactory::CreateRunManager(G4RunManagerType::Default);

  // 若提供第二个参数，则设置工作线程数
  if (argc == 3) {
    G4int nThreads = G4UIcommand::ConvertToInt(argv[2]);
    runManager->SetNumberOfThreads(nThreads);
  }

  // —— 注册强制用户初始化类（探测器、物理、动作）——
  DetectorConstruction* det = new DetectorConstruction;
  runManager->SetUserInitialization(det);

  PhysicsList* phys = new PhysicsList;
  runManager->SetUserInitialization(phys);

  runManager->SetUserInitialization(new ActionInitialization(/*det, phys*/));

  // —— 初始化可视化（仅在交互模式启用）——
  G4VisManager* visManager = nullptr;

  // 获取 UI 管理器指针
  G4UImanager* UImanager = G4UImanager::GetUIpointer();

  if (ui) {
    // 交互模式：启动可视化并进入交互会话
    visManager = new G4VisExecutive;
    visManager->Initialize();
    UImanager->ApplyCommand("/control/execute init_vis.mac");
    ui->SessionStart();
    delete ui;
  }
  else {
    // 批处理模式：先解析宏文件、把 source type / compartment 缓存到 RunAction
    // （不然 master 线程的 BeginOfRunAction 拿不到正确的源类型，会写出错的文件名）。
    ScanMacroForRunMeta(argv[1]);
    G4String command = "/control/execute ";
    G4String fileName = argv[1];
    UImanager->ApplyCommand(command + fileName);
  }

  // —— 作业结束：释放可视化与运行管理器 ——
  delete visManager;
  delete runManager;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
