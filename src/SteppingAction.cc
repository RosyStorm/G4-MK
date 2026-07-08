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

#include "DetectorConstruction.hh"
#include "EventAction.hh"
#include "G4ParticleDefinition.hh"
#include "G4RunManager.hh"
#include "G4Step.hh"
#include "G4StepPoint.hh"
#include "G4Track.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// 判断粒子是否在 α 链(He2+/He+/He0). Geant4-DNA 中 α 通过 G4DNAChargeDecrease 改
// 变电荷态时**会创建新粒子定义**(alpha→alpha+→alpha++→helium), 用 G4Alpha::Definition()
// 只匹配第一级, 漏掉后续链上粒子. 这里用粒子名判定覆盖完整 α 链.
namespace {
bool IsAlphaChain(const G4ParticleDefinition* pd) {
  if (!pd) return false;
  const G4String& name = pd->GetParticleName();
  // 覆盖 G4DNAChargeDecrease/Increase 链: alpha, alpha+, alpha++, helium
  return name == "alpha" || name == "alpha+" || name == "alpha++" || name == "helium";
}
}  // namespace

void SteppingAction::UserSteppingAction(const G4Step* step)
{
  // 任务6.2: 累加全局能量沉积(全部体积、全部粒子, 用于能量平衡校验)
  // 须在任何 kill 之前累加, 保证 killed 步的真实沉积也被计入
  if (fEventAction) fEventAction->AddEdep(step->GetTotalEnergyDeposit());

  G4Track* track = step->GetTrack();

  // ===== #1 加速: 出细胞(R_cell)且向外运动的粒子直接 kill =====
  // 核级/域级打分都在核内, 核外输运对结果无贡献。
  // 用 pos·dir>0 判定"向外", 以保留胞外源向内射入细胞的 α(交叉照射)。
  // 注意: 这会截断 α 射程, 故 α 射程验证(任务3.1)须把 /mygeom/killOutsideCell 设 false。
  const auto* det = static_cast<const DetectorConstruction*>(
    G4RunManager::GetRunManager()->GetUserDetectorConstruction());
  if (det && det->GetKillOutsideCell()) {
    // kill 半径: 默认 R_n(出核即杀, 更快且核打分无偏), 可切回 R_cell
    G4double Rkill =
      det->GetKillAtNucleus() ? det->GetNucleusRadius() : det->GetCellRadius();
    G4StepPoint* post = step->GetPostStepPoint();
    G4ThreeVector pos = post->GetPosition();
    if (pos.mag() > Rkill) {
      if (pos.dot(post->GetMomentumDirection()) > 0.) {
        track->SetTrackStatus(fStopAndKill);
        return;
      }
    }
  }

  // ===== P0 修复 #1: 出射边界步核内 edep 补加 =====
  // TrackerSD 默认按 volume 过滤: pre 在核内、post 在核外的 step, edep 归到 post (核外),
  // SD 不收到 hit, 导致 nucleusEdep (hit-sum) 漏掉这部分. 这里补加(保守按整 step 在核内).
  // 入射步 (pre 在核外 post 在核内) 由 SD 的 hit-sum 自然覆盖; 核内 step 同理.
  {
    const G4StepPoint* pre  = step->GetPreStepPoint();
    const G4StepPoint* post = step->GetPostStepPoint();
    const auto* preVol  = pre->GetTouchableHandle()->GetVolume();
    const auto* postVol = post->GetTouchableHandle()->GetVolume();
    if (preVol && postVol &&
        preVol->GetName() == "Nucleus" && postVol->GetName() != "Nucleus") {
      if (fEventAction)
        fEventAction->AddNucleusEdepBoundary(step->GetTotalEnergyDeposit());
    }
  }

  // ===== 任务3.1: α 链投影射程累计(射程验证时 kill 开关须关) =====
  if (!fEventAction) return;

  // 仅对 α 初级源记录射程：在事件首个初级粒子(应为 α)的第一步登记发射点
  if (!fEventAction->HaveVertex()) {
    if (track->GetParentID() == 0 && IsAlphaChain(track->GetParticleDefinition())) {
      fEventAction->SetPrimaryVertex(track->GetVertexPosition());
    }
    else {
      return;  // 非初级 α(如质子事件、δ 电子)在登记前不记录
    }
  }

  // 累计本事件 α 链相对发射点的最大位移 = 投影射程(直进，≈ CSDA 射程)
  // 注意: 只对 α 链上粒子(alpha/alpha+/alpha++/helium)计入. Geant4-DNA 中电荷交换
  //       会创建新粒子定义, 必须用粒子名判定覆盖全链, 否则漏掉 alpha+ 等后续粒子.
  //       δ 电子(G4Electron)的远端位置不被计入, 否则被电子拉宽射程直方图.
  if (!IsAlphaChain(track->GetParticleDefinition())) return;
  fEventAction->UpdateMaxRange(step->GetPostStepPoint()->GetPosition());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
