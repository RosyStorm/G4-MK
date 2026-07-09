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
/// \file SteppingAction.hh
/// \brief SteppingAction 类的定义：逐步动作
///
/// 任务3.1：记录初级 α 粒子停止时的 CSDA 射程(发射点→停止点直线距离)，
/// 用于与 NIST ASTAR 对比，验证 G4EmDNAPhysics_option2 的 α 输运。

#ifndef SteppingAction_h
#define SteppingAction_h 1

#include "G4UserSteppingAction.hh"

class EventAction;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

/// @brief 逐步动作类
///
/// 在每个 step 上执行三类处理：(1) 累加全局能量沉积（任务6.2 能量平衡）；
/// (2) 出细胞向外运动的粒子直接 kill（加速，核内打分不受影响）；
/// (3) 补加出射边界步的核内 edep（P0 修复 #1）；
/// (4) 累计 α 链(He2+→He+→He0)相对发射点的最大位移作为投影射程（任务3.1）。
class SteppingAction : public G4UserSteppingAction
{
  public:
    // ===== 构造与析构 =====
    /// 构造，绑定逐事件动作用于跨步累加。
    /// @param ea 关联的 EventAction 指针（可为空）
    explicit SteppingAction(EventAction* ea = nullptr) : fEventAction(ea) {}
    ~SteppingAction() override = default;  // 默认析构

    // ===== Geant4 强制重载接口 =====
    /// 每步回调：执行能量累加、kill 加速与射程累计。
    /// @param step 当前粒子步指针
    void UserSteppingAction(const G4Step*) override;

  private:
    EventAction* fEventAction = nullptr;  // 关联的逐事件动作（用于跨步累加）
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
