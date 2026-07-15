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
// * institutes,nor the agencies providing financial support to this  *
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
// *                                                                  *
// ********************************************************************
//
/// \file TrackingAction.hh
/// \brief TrackingAction 类的定义：径迹级用户动作
///
/// 任务X: 在每条 track 启动时立即调用 EventAction::EventParticleID 分类,
/// 防止 α 链衍生 track 在产生 hit 之前没被分类, 导致子粒子(δ 电子)无法
/// 继承 α 链根的 eid 而变成新 eid, 进而导致一个 α 被拆成 2-7 个"单事件
/// 粒子". 在 PreUserTrackingAction 钩子里强制分类, 解决此问题.

#ifndef TrackingAction_h
#define TrackingAction_h 1

#include "G4UserTrackingAction.hh"

class EventAction;

//....oooOO0OOOOooo........oooOO0OOOOooo........oooOO0OOOOooo........oooOO0OOOOooo......

/// @brief 径迹级用户动作类
///
/// 在 PreUserTrackingAction 中立即对每条新 track 调 EventParticleID,
/// 保证整条 track 链上每个节点都在 EventAction 的 fTrack2Event 缓存里,
/// 子粒子可以正确继承父粒子的 eid.
class TrackingAction : public G4UserTrackingAction
{
  public:
    // ===== 构造与析构 =====
    /// 构造, 绑定 EventAction.
    /// @param ea 关联的 EventAction 指针
    explicit TrackingAction(EventAction* ea = nullptr) : fEventAction(ea) {}
    ~TrackingAction() override = default;  // 默认析构

    // ===== Geant4 强制重载接口 =====
    /// 每条 track 启动前: 立即分类到 eid (强制写入缓存).
    void PreUserTrackingAction(const G4Track*) override;

  private:
    EventAction* fEventAction = nullptr;  // 关联的逐事件动作
};

#endif
