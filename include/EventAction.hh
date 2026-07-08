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
/// \brief Definition of the EventAction class
///
/// 任务3.1：逐事件累计 α 链(He2+→He+→He0 电荷交换子代)相对发射点的最大位移，
/// 即投影射程，填入 alphaRange 直方图。电荷交换会使原 α 径迹频繁"结束"，
/// 故必须跨整条链求最大位移，而非只跟初级。

#ifndef EventAction_h
#define EventAction_h 1

#include "G4ThreeVector.hh"
#include "G4UserEventAction.hh"
#include "globals.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class EventAction : public G4UserEventAction
{
  public:
    EventAction() = default;
    ~EventAction() override = default;

    void BeginOfEventAction(const G4Event*) override;
    void EndOfEventAction(const G4Event*) override;

    // 由 SteppingAction 调用
    void SetPrimaryVertex(const G4ThreeVector& v)
    {
      fPrimaryVertex = v;
      fHaveVertex = true;
    }
    G4bool HaveVertex() const { return fHaveVertex; }
    void UpdateMaxRange(const G4ThreeVector& pos);

    // 任务6.2: 全局能量沉积累加(SteppingAction 每步调用, 含全部体积/全部粒子)
    void AddEdep(G4double e) { fTotalEdep += e; }
    G4double GetTotalEdep() const { return fTotalEdep; }

    // P0 修复 #1: 出射边界步核内 edep 补加(SteppingAction 在 pre 在核内 post 在核外时调用)
    //   TrackerSD 默认按 volume 过滤, 跨边界出射 step 的 edep 归到 post volume(核外),
    //   SD 不收到 hit, 故 TrackerSD::EndOfEvent 中的 hit-sum nucleusEdep 会漏掉这部分.
    //   这里补加, 然后在 EndOfEvent 中把 hit-sum + boundary 求和得完整核内沉积.
    void AddNucleusEdepBoundary(G4double e) { fNucleusEdepBoundary += e; }
    G4double GetNucleusEdepBoundary() const { return fNucleusEdepBoundary; }

  private:
    G4ThreeVector fPrimaryVertex{};
    G4double fMaxRange = 0.;
    G4bool fHaveVertex = false;
    G4double fTotalEdep = 0.;  // 全局能量沉积(任务6.2 能量平衡)
    G4double fNucleusEdepBoundary = 0.;  // 边界步核内 edep 补加(P0 修复 #1)
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
