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
/// \file EventAction.hh
/// \brief EventAction 类的定义：逐事件动作
///
/// 任务3.1：逐事件累计 α 链(He2+→He+→He0 电荷交换子代)相对发射点的最大位移，
/// 即投影射程，填入 alphaRange 直方图。电荷交换会使原 α 径迹频繁"结束"，
/// 故必须跨整条链求最大位移，而非只跟初级。

#ifndef EventAction_h
#define EventAction_h 1

#include "G4ThreeVector.hh"
#include "G4UserEventAction.hh"
#include "globals.hh"

#include <map>

class G4Track;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

/// @brief 逐事件动作类
///
/// 在每个事件开始时重置累加器，事件结束时把 α 链投影射程填入直方图。
/// 同时提供全局能量沉积与核内边界步能量沉积的累加接口，供 SteppingAction 调用。
class EventAction : public G4UserEventAction
{
  public:
    // ===== 构造与析构 =====
    EventAction() = default;            // 默认构造
    ~EventAction() override = default;  // 默认析构

    // ===== Geant4 强制重载接口 =====
    void BeginOfEventAction(const G4Event*) override;  // 事件开始：重置累加器
    void EndOfEventAction(const G4Event*) override;     // 事件结束：填充射程直方图

    // ===== 由 SteppingAction 调用的接口 =====
    /// 登记初级 α 发射点（由事件首个初级粒子的第一步调用）。
    /// @param v 初级粒子的产生位置（vertex position）
    void SetPrimaryVertex(const G4ThreeVector& v)
    {
      fPrimaryVertex = v;
      fHaveVertex = true;
    }
    /// 是否已登记发射点。
    /// @return 已登记返回 true
    G4bool HaveVertex() const { return fHaveVertex; }

    /// 更新 α 链相对发射点的最大位移（投影射程）。
    /// @param pos 当前 α 链粒子的位置
    void UpdateMaxRange(const G4ThreeVector& pos);

    // ===== 全局能量沉积累加（任务6.2：能量平衡）=====
    // 含全部体积、全部粒子，由 SteppingAction 每步调用
    /// 累加全局能量沉积。
    /// @param e 本步总能量沉积
    void AddEdep(G4double e) { fTotalEdep += e; }
    /// 取本事件全局能量沉积累加值。
    /// @return 全局能量沉积（内部单位）
    G4double GetTotalEdep() const { return fTotalEdep; }

    // ===== 出射边界步核内 edep 补加（P0 修复 #1）=====
    // TrackerSD 默认按 volume 过滤，跨边界出射 step 的 edep 归到 post volume(核外),
    // SD 不收到 hit，故 TrackerSD::EndOfEvent 中的 hit-sum nucleusEdep 会漏掉这部分。
    // 这里补加，然后在 EndOfEvent 中把 hit-sum + boundary 求和得完整核内沉积。
    /// 补加出射边界步的核内能量沉积。
    /// @param e 该出射边界步的总能量沉积
    void AddNucleusEdepBoundary(G4double e) { fNucleusEdepBoundary += e; }
    /// 取边界步核内能量沉积补加值。
    /// @return 边界步核内能量沉积（内部单位）
    G4double GetNucleusEdepBoundary() const { return fNucleusEdepBoundary; }

    // ===== 按粒子分组(单事件 f_{n,1}/f_{d,1}; 路线2 用, 路线1 退化为 1 组)=====
    /// 取该 track 所属"单事件粒子"的 ID。
    /// 分类: 创建过程为 fDecay(RadioactiveDecay 产物: α/β/反冲) 或无创建者(初级) → 新 eventID;
    ///       电离次级(δ 电子) → 继承母粒子的 eventID。结果缓存于 fTrack2Event。
    /// @param track 当前径迹
    /// @return 单事件粒子 ID(≥0)
    G4int EventParticleID(const G4Track* track) const;
    /// 取某 eventID 对应粒子的 PDG 编码(填 ntuple 用)。
    /// @param eid 单事件粒子 ID
    /// @return PDG 编码(未登记返回 0)
    G4int EventParticlePDG(G4int eid) const
    {
      auto it = fEvent2PDG.find(eid);
      return it != fEvent2PDG.end() ? it->second : 0;
    }
    /// 取某 eventID 对应粒子的动能(填 ntuple alphaE_MeV 用)。
    /// @param eid 单事件粒子 ID
    /// @return 产生时动能(未登记返回 0)
    G4double EventParticleKE(G4int eid) const
    {
      auto it = fEvent2KE.find(eid);
      return it != fEvent2KE.end() ? it->second : 0.;
    }

  private:

  private:
    // ===== 逐事件累加量 =====
    G4ThreeVector fPrimaryVertex{};              // 初级 α 发射点
    G4double fMaxRange = 0.;                     // α 链相对发射点的最大位移（投影射程）
    G4bool fHaveVertex = false;                  // 是否已登记发射点
    G4double fTotalEdep = 0.;                    // 全局能量沉积（任务6.2 能量平衡）
    G4double fNucleusEdepBoundary = 0.;          // 边界步核内 edep 补加（P0 修复 #1）

    // ===== 按粒子分组(单事件; 路线2)=====
    mutable std::map<G4int, G4int> fTrack2Event;  // trackID → 单事件粒子 ID
    mutable std::map<G4int, G4int> fEvent2PDG;    // 单事件粒子 ID → PDG 编码
    mutable std::map<G4int, G4double> fEvent2KE;   // 单事件粒子 ID → 产生时动能
    mutable G4int fNextEventID = 0;                // 下一个单事件粒子 ID
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
