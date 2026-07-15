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
/// \file EventAction.cc
/// \brief EventAction 类的实现：事件级用户动作
///
/// 在事件开始时重置事件级累积量，事件进行中累积能量沉积与 α 投影射程，
/// 事件末将结果（如 α 射程直方图）填入分析管理器。

#include "EventAction.hh"

#include "G4AnalysisManager.hh"
#include "G4Event.hh"
#include "G4ParticleDefinition.hh"
#include "G4ProcessType.hh"
#include "G4SystemOfUnits.hh"
#include "G4Track.hh"
#include "G4VProcess.hh"
#include "globals.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::BeginOfEventAction(const G4Event* /*anEvent*/)
{
  /// 事件开始时重置所有事件级累积量。
  /// @param anEvent 当前事件（未使用）

  // —— 重置事件级累积量 ——
  fMaxRange = 0.;                       // α 最大位移（投影射程）清零
  fHaveVertex = false;                  // 尚未登记发射点
  fPrimaryVertex = G4ThreeVector();     // 初级 α 发射点清零
  fTotalEdep = 0.;                      // 任务6.2: 全局能量沉积清零
  fNucleusEdepBoundary = 0.;            // P0 修复 #1: 边界步核内 edep 累加器清零
  // 按粒子分组(单事件): 清映射表
  fTrack2Event.clear();
  fEvent2PDG.clear();
  fEvent2KE.clear();
  fNextEventID = 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4int EventAction::EventParticleID(const G4Track* track) const
{
  /// 取该 track 所属"单事件粒子"的 ID(带缓存)。
  ///
  /// 分类规则(任务X更新):
  ///   1. 初级(枪直接产生, creator==nullptr)         → 新单事件粒子
  ///   2. 真正的放射性衰变产物(creator->GetProcessName()=="RadioactiveDecay")
  ///                                               → 新单事件粒子
  ///   3. 其他所有粒子(包括 α 链衍生 α⁺/α²⁺/He, G4-DNA 电荷交换
  ///      alpha_G4DNAChargeDecrease 等)             → 继承母粒子 ID
  ///
  /// 任务X设计动机: 1 个 5.83 MeV α 在核内输运时, G4-DNA 通过电荷交换
  ///   过程(alpha_G4DNAChargeDecrease)产生 alpha+ track, 再经
  ///   alpha+_G4DNAChargeIncrease 回到 alpha. 这些衍生 track 都属于
  ///   "同一个 α 在核内的同一次穿越", 应整体视为 1 个单事件粒子
  ///   (Kellerer 1975 microdosimetry: 一次能量沉积事件 = 1 个带电粒子
  ///   在 site 内的完整穿越). 旧代码因 ProcessType 误判把每个衍生
  ///   track 当作"新单事件粒子", 导致每事件 1-7 个 α 行.
  ///
  /// @param track 当前径迹
  /// @return 单事件粒子 ID(≥0)

  if (!track) return -1;
  G4int tid = track->GetTrackID();
  auto it = fTrack2Event.find(tid);
  if (it != fTrack2Event.end()) return it->second;  // 已缓存

  G4int eid = -1;
  const G4VProcess* creator = track->GetCreatorProcess();
  G4bool isPrimary = (creator == nullptr);
  G4bool isRealRadioactiveDecay = (!isPrimary && creator &&
                                   creator->GetProcessName() == "RadioactiveDecay");

  if (isPrimary || isRealRadioactiveDecay) {
    // —— 新单事件粒子 ——
    eid = fNextEventID++;
    if (track->GetParticleDefinition()) {
      fEvent2PDG[eid] = track->GetParticleDefinition()->GetPDGEncoding();
    }
    fEvent2KE[eid] = track->GetVertexKineticEnergy();
  }
  else {
    // —— 任何次级(包括 α 链衍生)继承母粒子 ID ——
    G4int pid = track->GetParentID();
    auto pit = fTrack2Event.find(pid);
    if (pit != fTrack2Event.end()) {
      eid = pit->second;
    }
    else {
      // 母粒子未被分类(未在核内 hit). 新建一个 eventID 但记为继承 ——
      // 这里不能简单地 fNextEventID++, 否则 e⁻ 抢占 eid 顺序会破坏
      // α 链继承. 正确做法: 沿父链向上追到已分类节点.
      G4int curPid = pid;
      G4int safety = 0;
      while (curPid > 0 && safety++ < 64) {
        auto cit = fTrack2Event.find(curPid);
        if (cit != fTrack2Event.end()) {
          eid = cit->second;
          break;
        }
        // 父粒子未被分类: 拿不到父 PDG, 此处无法判断继续向上 ——
        // G4 没有暴露 track-by-parentId 查 PDG 的接口, 所以 fallback
        // 到新建 eid (虽然可能与 α 链共享 eid 不一致, 但仅发生在
        // 母粒子完全没 hit 时, 概率极低).
        break;
      }
      if (eid < 0) {
        eid = fNextEventID++;
        if (track->GetParticleDefinition()) {
          fEvent2PDG[eid] = track->GetParticleDefinition()->GetPDGEncoding();
        }
        fEvent2KE[eid] = track->GetVertexKineticEnergy();
      }
    }
  }
  fTrack2Event[tid] = eid;
  return eid;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::UpdateMaxRange(const G4ThreeVector& pos)
{
  /// 更新 α 相对发射点的最大位移（投影射程）。
  /// @param pos 当前步后位置
  /// @return 无（结果写入成员 fMaxRange）

  // 尚未登记发射点：直接返回，无需更新
  if (!fHaveVertex) return;

  // —— 计算当前位置相对发射点的距离 ——
  G4double distance = (pos - fPrimaryVertex).mag();

  // —— 保留最大值（即投影射程）——
  if (distance > fMaxRange) fMaxRange = distance;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::EndOfEventAction(const G4Event* /*anEvent*/)
{
  /// 事件末动作：若有有效发射点，将 α 投影射程填入直方图。
  /// @param anEvent 当前事件（未使用）

  // —— 检查是否有有效发射点与射程 ——
  if (fHaveVertex && fMaxRange > 0.) {
    // —— 填充 α 射程直方图 ——
    auto* analysisManager = G4AnalysisManager::Instance();
    G4int histogramId = analysisManager->GetH1Id("alphaRange");
    if (histogramId >= 0) {
      analysisManager->FillH1(histogramId, fMaxRange / um);
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
