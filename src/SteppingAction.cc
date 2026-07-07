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
/// \file SteppingAction.cc
/// \brief Implementation of the SteppingAction class
///
/// 任务3.1：累计 α 链(He2+ 经 G4DNAChargeDecrease 变 He+/He0)相对发射点的
/// 最大位移，作为投影射程填入 alphaRange。电荷交换频繁(平均自由程~um)，
/// 故必须跨整条链求最大位移，而非只跟初级径迹。

#include "SteppingAction.hh"

#include "EventAction.hh"
#include "G4Alpha.hh"
#include "G4Step.hh"
#include "G4StepPoint.hh"
#include "G4Track.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SteppingAction::UserSteppingAction(const G4Step* step)
{
  if (!fEventAction) return;
  const G4Track* track = step->GetTrack();

  // 仅对 α 初级源记录射程：在事件首个初级粒子(应为 α)的第一步登记发射点
  if (!fEventAction->HaveVertex()) {
    if (track->GetParentID() == 0 &&
        track->GetParticleDefinition() == G4Alpha::Definition()) {
      fEventAction->SetPrimaryVertex(track->GetVertexPosition());
    }
    else {
      return;  // 非初级 α(如质子事件、δ 电子)在登记前不记录
    }
  }

  // 累计本事件 α 链相对发射点的最大位移 = 投影射程(直进，≈ CSDA 射程)
  fEventAction->UpdateMaxRange(step->GetPostStepPoint()->GetPosition());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
