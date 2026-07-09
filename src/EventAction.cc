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
#include "G4SystemOfUnits.hh"
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
