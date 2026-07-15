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
/// \file TrackingAction.cc
/// \brief TrackingAction 类的实现：径迹级用户动作
///
/// 任务X: 在每条 track 启动时立即调 EventAction::EventParticleID 强制分类.

#include "TrackingAction.hh"

#include "EventAction.hh"
#include "G4Event.hh"
#include "G4RunManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4Track.hh"

//....oooOO0OOOOooo........oooOO0OOOOooo........oooOO0OOOOooo........oooOO0OOOOooo......

void TrackingAction::PreUserTrackingAction(const G4Track* track)
{
  /// 每条 track 启动前: 立即按 EventAction::EventParticleID 规则分类,
  /// 强制写入 fTrack2Event 缓存. 这样 α 链衍生 track 即使没在核内 hit,
  /// 也会被先分类, 后续 δ 电子可以正确继承 α 链根 eid.

  if (fEventAction && track) {
    fEventAction->EventParticleID(track);

    // 诊断已移除(确认: 多 eid 是 Hf-177m 内部转换电子, 是真实物理)
    (void)0;
  }
}
