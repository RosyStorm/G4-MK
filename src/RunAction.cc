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
/// \file RunAction.cc
/// \brief RunAction 类的实现：Run 级动作
///
/// 负责初始化 Geant4 分析管理器(ROOT 输出)、创建直方图与 ntuple、
/// 在 Run 结束时打印微剂量学统计（频率均值/剂量均值等）并写入 ROOT 文件。

#include "RunAction.hh"

#include "G4AnalysisManager.hh"
#include "G4RunManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4UnitsTable.hh"
#include "globals.hh"

#include "PrimaryGeneratorAction.hh"  // 取源类型/区室 → 动态文件名

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// —— 静态缓存：源类型 / 区室，由 microtrack.cc 解析宏文件后注入 ——
namespace {
G4String gRunSourceType   = "ac225_decay";  // 兜底
G4String gRunCompartment  = "Membrane";     // 兜底
}

void RunAction::SetRunMeta(G4String sourceType, G4String compartment)
{
  /// 由 microtrack.cc 在 ApplyCommand 之前调用，确保 master 线程也能拿到正确的
  /// 源类型与区室（workder 线程 PGA 不在 master 上注册）。
  gRunSourceType  = sourceType;
  gRunCompartment = compartment;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::RunAction()
{
  /// 构造函数：创建分析管理器并定义全部直方图与 ntuple 的分箱方案。

  // —— 创建分析管理器 ——
  auto analysisManager = G4AnalysisManager::Instance();

  analysisManager->SetVerboseLevel(1);
  analysisManager->SetNtupleMerging(true);

  // ===== 定义各类直方图的分箱方案 =====

  // —— 动能入射/出射直方图的线性分箱 ——
  const G4double kinEmin = 0.;            // 动能下限
  const G4double kinEmax = 300.;          // 动能上限
  const G4double kinEbinWidth = 0.1;      // bin 宽
  const G4int kinEbins =
    static_cast<G4int>((kinEmax - kinEmin) / kinEbinWidth);

  // —— 能量直方图的对数分箱 ——
  const G4int minLog10E = -4.;            // log10 下限
  const G4int maxLog10E = +4.;            // log10 上限
  const G4int nBinsLog10E = (maxLog10E - minLog10E) * 50;  // bin 数(每量级 50 bin)

  G4double binsLog10E[nBinsLog10E + 1];
  G4double binWidthLog10E =
    static_cast<G4double>((maxLog10E - minLog10E)) / nBinsLog10E;

  for (G4int ii = 0; ii <= nBinsLog10E; ii++) {
    binsLog10E[ii] = std::pow(10., binWidthLog10E * ii + minLog10E);
  }

  std::vector<G4double> vecBinsLog10E(binsLog10E, binsLog10E + nBinsLog10E + 1);

  // —— 加权计数直方图的对数分箱 ——
  const G4int minLog10W = 0;
  const G4int maxLog10W = +6;
  const G4int nBinsLog10W = (maxLog10W - minLog10W) * 50;

  G4double binsLog10W[nBinsLog10W + 1];
  G4double binWidthLog10W =
    static_cast<G4double>((maxLog10W - minLog10W)) / nBinsLog10W;

  for (G4int ii = 0; ii <= nBinsLog10W; ii++) {
    binsLog10W[ii] = std::pow(10., binWidthLog10W * ii + minLog10W);
  }

  std::vector<G4double> vecBinsLog10W(binsLog10W, binsLog10W + nBinsLog10W + 1);

  // —— 计数直方图的线性分箱 ——
  // 注: α 等高 LET 粒子在核内可产生数万~数十万 hit, 上限需足够大以免溢出致 mean 失真
  G4int minCount = 0;
  G4int maxCount = 500000;
  G4int nBinsCount = (maxCount - minCount) / 100;

  // ===== 创建直方图 =====

  // —— 能量谱直方图(对数分箱) ——
  analysisManager->CreateH1(
    "fe",
    "Energy imparted per event [keV] (log binning)",
    vecBinsLog10E);

  analysisManager->CreateH1(
    "efe",
    "Weighted energy imparted per event [keV] (log binning)",
    vecBinsLog10E);

  analysisManager->CreateH1(
    "e2fe",
    "Squared-weighted energy imparted per event [keV] (log binning)",
    vecBinsLog10E);

  // —— 线能直方图(对数分箱) ——
  analysisManager->CreateH1("fy", "Lineal energy [keV/um] (log binning)",
                            vecBinsLog10E);

  analysisManager->CreateH1(
    "yfy",
    "Dose-weighted lineal energy [keV/um] (log binning)", vecBinsLog10E);

  analysisManager->CreateH1(
    "y2fy",
    "Squared-weighted lineal energy [keV/um] (log binning)", vecBinsLog10E);

  // —— 比能直方图(对数分箱) ——
  analysisManager->CreateH1(
    "fz",
    "Single-event specific energy [Gy] (log binning)",
    vecBinsLog10E);

  analysisManager->CreateH1(
    "zfz",
    "Dose-weighted single-event specific energy [Gy] (log binning)",
    vecBinsLog10E);

  analysisManager->CreateH1(
    "z2fz",
    "Squared-weighted single-event specific energy [Gy] (log binning)",
    vecBinsLog10E);

  // —— 命中计数直方图(线性分箱) ——
  analysisManager->CreateH1(
    "Nsel", "Number of selectable hits per event",
    nBinsCount, minCount, maxCount);

  analysisManager->CreateH1(
    "Nsite", "Number of hits in site", nBinsCount,
    minCount, maxCount);

  analysisManager->CreateH1("Nint",
    "Number of selectable hits in site", nBinsCount,
    minCount, maxCount);

  // —— 核边界入射/出射动能直方图 ——
  analysisManager->CreateH1("KinE_in", "Kinetic energy at the entrance [MeV]",
                            kinEbins, kinEmin, kinEmax);
  analysisManager->CreateH1("KinE_out", "Kinetic energy at the exit [MeV]",
                            kinEbins, kinEmin, kinEmax);

  // —— 域内命中数 vs 加权能量沉积 (2D) ——
  analysisManager->CreateH2(
    "Nsite_vs_e",
    "Number of hits in site vs energy imparted [keV] (log-log)",
    vecBinsLog10E, vecBinsLog10W);

  // 任务3.1：初级 α 的 CSDA 射程(与 NIST ASTAR 对比)
  analysisManager->CreateH1(
    "alphaRange", "Primary alpha CSDA range [um] (task 3.1)",
    240, 0., 120.);

  // 任务4.1：核单事件比能 z_n 谱(与 z_d 同对数分箱)，→ z̄_{n,F}, z̄_{n,D} 喂 SMK 式(24)
  analysisManager->CreateH1(
    "fzn", "Single-event nucleus specific energy z_n [Gy] (log binning)",
    vecBinsLog10E);
  analysisManager->CreateH1(
    "znfzn", "Dose-weighted z_n [Gy] (log binning)", vecBinsLog10E);
  analysisManager->CreateH1(
    "z2nfzn", "Squared-weighted z_n [Gy] (log binning)", vecBinsLog10E);

  // —— 任务4.2：逐事件配对 ntuple —— 每事件一行，hit/miss 都填 ——
  // 后处理(ROOT)可算 z_d/z_n 边缘谱、联合分布、条件 f(z_d|z_n)、多事件卷积与 SMK 存活曲线
  analysisManager->CreateNtuple("events", "Per-event microdosimetry (task 4.2)");
  analysisManager->CreateNtupleIColumn("eventID");        // 0
  analysisManager->CreateNtupleDColumn("alphaE_MeV");     // 1 初级 α 动能
  analysisManager->CreateNtupleDColumn("edep_d_keV");     // 2 域内能量沉积 ε_d
  analysisManager->CreateNtupleDColumn("z_d_Gy");         // 3 域比能 z_d
  analysisManager->CreateNtupleDColumn("edep_n_keV");     // 4 核内总能量沉积 ε_n
  analysisManager->CreateNtupleDColumn("z_n_Gy");         // 5 核比能 z_n
  analysisManager->CreateNtupleDColumn("weight");         // 6 域抽样权 w=Nsel/Nint
  analysisManager->CreateNtupleIColumn("nHsel");          // 7
  analysisManager->CreateNtupleIColumn("nHsite");         // 8
  analysisManager->CreateNtupleIColumn("nHint");          // 9
  analysisManager->CreateNtupleIColumn("hitFlag");        // 10  1=命中核, 0=miss
  analysisManager->CreateNtupleIColumn("compartment");    // 11  0=Nuc,1=Cyt,2=Mem,3=Ext
  analysisManager->CreateNtupleDColumn("edep_total_keV"); // 12 全局能量沉积(任务6.2)
  analysisManager->FinishNtuple();

  // 按粒子分组 ntuple (id=1, 每个单事件粒子一行) —— 路线2 出 f_{n,1}/f_{d,1}
  // 列 0-12 与 events (id=0) 完全一致(共用分析代码); 13-14 为路线2 额外列
  analysisManager->CreateNtuple("single_events",
    "Per-particle single-event (task 7.1: route-2 f_{n,1}/f_{d,1})");
  analysisManager->CreateNtupleIColumn("eventID");         // 0  (同 events)
  analysisManager->CreateNtupleDColumn("alphaE_MeV");      // 1  该粒子产生时动能
  analysisManager->CreateNtupleDColumn("edep_d_keV");      // 2
  analysisManager->CreateNtupleDColumn("z_d_Gy");          // 3
  analysisManager->CreateNtupleDColumn("edep_n_keV");      // 4
  analysisManager->CreateNtupleDColumn("z_n_Gy");          // 5
  analysisManager->CreateNtupleDColumn("weight");          // 6
  analysisManager->CreateNtupleIColumn("nHsel");           // 7  该粒子核内 hit 数
  analysisManager->CreateNtupleIColumn("nHsite");           // 8  该粒子域内 hit 数
  analysisManager->CreateNtupleIColumn("nHint");            // 9
  analysisManager->CreateNtupleIColumn("hitFlag");          // 10 恒为 1(有沉积才入此表)
  analysisManager->CreateNtupleIColumn("compartment");      // 11
  analysisManager->CreateNtupleDColumn("edep_total_keV");  // 12 该粒子核沉积(≈edep_n)
  analysisManager->CreateNtupleIColumn("eventParticleID");  // 13 [路线2额外] 单事件粒子 ID
  analysisManager->CreateNtupleIColumn("pdg");              // 14 [路线2额外] PDG 编码
  analysisManager->FinishNtuple();

  // 任务X: 仅记录由 ac225_single_decay 模式下的 He-4 α 粒子单事件
  // 与 single_events 完全相同的列结构, 但仅在 ac225_single_decay 模式且 PDG == 1000020040
  // (He-4 α) 时由 TrackerSD::EndOfEvent 条件性写入
  analysisManager->CreateNtuple("alpha_events",
    "Per-He4-alpha single-event (ac225_single_decay only)");
  analysisManager->CreateNtupleIColumn("eventID");         // 0
  analysisManager->CreateNtupleDColumn("alphaE_MeV");      // 1  α 产生时动能
  analysisManager->CreateNtupleDColumn("edep_d_keV");      // 2
  analysisManager->CreateNtupleDColumn("z_d_Gy");          // 3
  analysisManager->CreateNtupleDColumn("edep_n_keV");      // 4
  analysisManager->CreateNtupleDColumn("z_n_Gy");          // 5
  analysisManager->CreateNtupleDColumn("weight");          // 6
  analysisManager->CreateNtupleIColumn("nHsel");           // 7
  analysisManager->CreateNtupleIColumn("nHsite");          // 8
  analysisManager->CreateNtupleIColumn("nHint");           // 9
  analysisManager->CreateNtupleIColumn("hitFlag");         // 10
  analysisManager->CreateNtupleIColumn("compartment");     // 11
  analysisManager->CreateNtupleDColumn("edep_total_keV");  // 12
  analysisManager->CreateNtupleIColumn("eventParticleID"); // 13
  analysisManager->CreateNtupleIColumn("pdg");             // 14
  analysisManager->FinishNtuple();

  // 任务Y: beta_events (id=3)
  // 与 alpha_events 完全相同的列结构, 但仅在 lu177_decay 模式且 PDG == 11
  // (e⁻ β⁻) 时由 TrackerSD::EndOfEvent 条件性写入。列 1 改名为 betaE_MeV
  // 以反映该 ntuple 仅含 β⁻ 物理量(而 alpha_events 是 α 物理量)。
  // 其它 ntuple(events/single_events/alpha_events)保留 alphaE_MeV 列名以避免
  // 改动既有分析脚本。
  analysisManager->CreateNtuple("beta_events",
    "Per-beta single-event (lu177_decay only, PDG==11)");
  analysisManager->CreateNtupleIColumn("eventID");         // 0
  analysisManager->CreateNtupleDColumn("betaE_MeV");       // 1  β⁻ 产生时动能
  analysisManager->CreateNtupleDColumn("edep_d_keV");      // 2
  analysisManager->CreateNtupleDColumn("z_d_Gy");          // 3
  analysisManager->CreateNtupleDColumn("edep_n_keV");      // 4
  analysisManager->CreateNtupleDColumn("z_n_Gy");          // 5
  analysisManager->CreateNtupleDColumn("weight");          // 6
  analysisManager->CreateNtupleIColumn("nHsel");           // 7
  analysisManager->CreateNtupleIColumn("nHsite");          // 8
  analysisManager->CreateNtupleIColumn("nHint");           // 9
  analysisManager->CreateNtupleIColumn("hitFlag");         // 10
  analysisManager->CreateNtupleIColumn("compartment");     // 11
  analysisManager->CreateNtupleDColumn("edep_total_keV");  // 12
  analysisManager->CreateNtupleIColumn("eventParticleID"); // 13
  analysisManager->CreateNtupleIColumn("pdg");             // 14
  analysisManager->FinishNtuple();
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::~RunAction() = default;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::BeginOfRunAction(const G4Run* /*aRun*/)
{
  /// Run 开始时打开 ROOT 输出文件。
  /// @param aRun 当前 Run 对象（未使用）

  // 通知 runManager 保存随机数种子（当前注释掉）
  // G4RunManager::GetRunManager()->SetRandomNumberStore(true);

  // —— 获取分析管理器 ——
  auto analysisManager = G4AnalysisManager::Instance();

  // —— 根据源类型 + 区室动态生成文件名 ——
  // 统一放在 data/ 下, 按 {源类型}_{区室}.root 命名, 避免不同配置互相覆盖
  // 注：master 线程上 GetUserPrimaryGeneratorAction() 返回 nullptr（MT 模式下
  //     PGA 由 ActionInitialization::Build() 在 worker 线程里创建, master 只有 RunAction）。
  //     因此这里直接用 microtrack.cc 从 argv[1] 解析后写入的静态缓存 gRunSourceType /
  //     gRunCompartment，避免 worker 上 PGA 与 master 上文件名不一致。设兜底与默认源类型一致。
  G4String sourceType  = gRunSourceType;
  G4String compartment = gRunCompartment;

  G4String fileName;
  if (sourceType == "ac225_decay") {
    fileName = "data/ac225_phy_decay_" + compartment + ".root";
  }
  else if (sourceType == "ac225_single_decay") {
    fileName = "data/ac225_single_decay_" + compartment + ".root";
  }
  else if (sourceType == "lu177_decay") {
    fileName = "data/lu177_phy_decay_" + compartment + ".root";
  }
  else if (sourceType == "ac225") {
    fileName = "data/4alpha_" + compartment + ".root";
  }
  else if (sourceType == "proton") {
    fileName = "data/proton_" + compartment + ".root";
  }
  else if (sourceType == "alpha") {
    fileName = "data/alpha.root";   // 单 α 验证, 无区室
  }
  else if (sourceType == "carbon") {
    fileName = "data/carbon_" + compartment + ".root";  // 碳离子验证, 按区室命名
  }
  else if (sourceType == "am241_decay") {
    fileName = "data/am241_phy_decay_" + compartment + ".root";  // Am-241 完整衰变链, 按区室命名
  }
  else {
    fileName = "data/microtrack.root";  // 兜底
  }

  analysisManager->OpenFile(fileName);
  G4cout << "输出文件: " << fileName << G4endl;
  G4cout << "使用 " << analysisManager->GetType() << " 分析管理器" << G4endl;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::EndOfRunAction(const G4Run* /*aRun*/)
{
  /// Run 结束时打印直方图统计并写入/关闭 ROOT 文件。
  /// @param aRun 当前 Run 对象（未使用）

  auto analysisManager = G4AnalysisManager::Instance();

  // —— 仅主线程打印整 Run 统计（且直方图 1 存在时）——
  if (IsMaster() && analysisManager->GetH1(1)) {

    G4cout << G4endl;

    G4cout << "----> 整个 Run 的直方图统计: "
           << G4endl << G4endl;

    // —— 打印直方图统计 ——
    G4cout << "  单事件授与能:\n"
           << "  ----------------------------"
           << G4endl << G4endl;

    G4cout << "    频率均值: \\varepsilon_{1,F} = "
           << analysisManager->GetH1(0)->mean() << " keV "
           << " (rms = " << analysisManager->GetH1(0)->rms() << " keV)"
           << G4endl;

    G4cout << "    剂量均值: \\varepsilon_{1,D} = "
           << analysisManager->GetH1(1)->mean() << " keV "
           << " (rms = " << analysisManager->GetH1(1)->rms() << " keV)"
           << G4endl;

    G4cout << "    加权平方直方图: 均值 = "
           << analysisManager->GetH1(2)->mean() << " keV "
           << " (rms = " << analysisManager->GetH1(2)->rms() << " keV)"
           << G4endl;

    G4cout << G4endl
           << "  线能:\n"
           << "  -------------"
           << G4endl << G4endl;

    G4cout << "    频率均值: y_F = "
           << analysisManager->GetH1(3)->mean() << " keV/um "
           << " (rms = " << analysisManager->GetH1(3)->rms() << " keV/um)"
           << G4endl;

    G4cout << "    剂量均值: y_D = "
           << analysisManager->GetH1(4)->mean() << " keV/um "
           << " (rms = " << analysisManager->GetH1(4)->rms() << " keV/um)"
           << G4endl;

    G4cout << "    加权平方直方图: 均值 = "
           << analysisManager->GetH1(5)->mean() << " keV/um "
           << " (rms = " << analysisManager->GetH1(5)->rms() << " keV/um)"
           << G4endl;

    G4cout << G4endl
           << "  单事件比能:\n"
           << "  ----------------------------"
           << G4endl << G4endl;

    G4cout << "    频率均值: z_{1,F} = "
           << analysisManager->GetH1(6)->mean() << " Gy "
           << " (rms = " << analysisManager->GetH1(6)->rms() << " Gy)"
           << G4endl;

    G4cout << "    剂量均值: z_{1,D} = "
           << analysisManager->GetH1(7)->mean() << " Gy "
           << " (rms = " << analysisManager->GetH1(7)->rms() << " Gy)"
           << G4endl;

    G4cout << "    加权平方直方图: 均值 = "
           << analysisManager->GetH1(8)->mean() << " Gy"
           << " (rms = " << analysisManager->GetH1(8)->rms() << " Gy)"
           << G4endl;

    // 任务4.1：核级 z_n (→ z̄_{n,D} 喂 SMK)
    G4int idZnF = analysisManager->GetH1Id("fzn");
    G4int idZnD = analysisManager->GetH1Id("znfzn");
    if (idZnF >= 0 && analysisManager->GetH1(idZnF) && idZnD >= 0
        && analysisManager->GetH1(idZnD)) {
      G4cout << G4endl
             << "  核单事件比能 (z_n, 任务4.1):\n"
             << "  ---------------------------------------------------"
             << G4endl << G4endl;
      G4cout << "    频率均值: z_{n,F} = "
             << analysisManager->GetH1(idZnF)->mean() << " Gy "
             << " (rms = " << analysisManager->GetH1(idZnF)->rms() << " Gy)"
             << G4endl;
      G4cout << "    剂量均值:      z_{n,D} = "
             << analysisManager->GetH1(idZnD)->mean() << " Gy "
             << " (rms = " << analysisManager->GetH1(idZnD)->rms() << " Gy)"
             << G4endl;
    }

    G4cout << G4endl
           << "  每事件命中数:\n"
           << "  ------------------------"
           << G4endl << G4endl;

    G4cout << "    可用于域随机放置的, N_{sel}: 均值 = "
           << analysisManager->GetH1(9)->mean()
           << " (rms = " << analysisManager->GetH1(9)->rms() << ")"
           << G4endl;

    G4cout << "    域内的, N_{site}: 均值 = "
           << analysisManager->GetH1(10)->mean()
           << " (rms = " << analysisManager->GetH1(10)->rms() << ")"
           << G4endl;

    G4cout << "    域内且可用于域随机放置的,"
           << " N_{int}: 均值 = " << analysisManager->GetH1(11)->mean()
           << " (rms = " << analysisManager->GetH1(11)->rms() << ")"
           << G4endl;

    G4cout << G4endl
           << "  初级粒子动能:\n"
           << "  --------------------------------------"
           << G4endl << G4endl;

    G4cout << "    入射核时, T_{in}: 均值 = "
           << G4BestUnit(analysisManager->GetH1(12)->mean(), "Energy")
           << " (rms = "
           << G4BestUnit(analysisManager->GetH1(12)->rms(), "Energy") << ")"
           << G4endl;

    G4cout << "    出射核时, T_{out}: 均值 = "
           << G4BestUnit(analysisManager->GetH1(13)->mean(), "Energy")
           << " (rms = "
           << G4BestUnit(analysisManager->GetH1(13)->rms(), "Energy") << ")"
           << G4endl;

    G4cout << G4endl;
  }

  // —— 保存直方图与 ntuple 到文件 ——
  analysisManager->Write();
  analysisManager->CloseFile();

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
