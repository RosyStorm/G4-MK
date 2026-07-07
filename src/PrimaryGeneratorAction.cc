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
/// \file PrimaryGeneratorAction.cc
/// \brief Implementation of the PrimaryGeneratorAction class
///
/// 支持两种源类型(/source/type)：
///   - proton : 原质子枪(基线对比)，位置/方向由成员决定，能量由 /gun/energy 给
///   - ac225  : Ac-225 α 源，从指定区室(/source/compartment)均匀随机点各向同性
///              发射一个 α，能量按 Ac-225 衰变链抽样(一个 α = 一个事件，契合 MK 单事件定义)
///   区室: Nucleus(核内) | Cytoplasm(质内) | Membrane(膜面,默认) | Extracellular(胞外)

#include "PrimaryGeneratorAction.hh"

#include "G4IonTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4RunManager.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"
#include "globals.hh"

#include <cmath>

#include "DetectorConstruction.hh"
#include "PrimaryGeneratorMessenger.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction::PrimaryGeneratorAction()
  : G4VUserPrimaryGeneratorAction()
{
  // Create a messenger for this class
  fGunMessenger = std::make_unique<PrimaryGeneratorMessenger>(this);

  // Create a particle gun
  G4int nOfParticles = 1;
  fParticleGun = std::make_unique<G4ParticleGun>(nOfParticles);

  // Set default particles.
  fParticle = G4ParticleTable::GetParticleTable()->FindParticle("proton");
  fAlpha = G4ParticleTable::GetParticleTable()->FindParticle("alpha");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction::~PrimaryGeneratorAction() = default;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  if (fSourceType == "ac225") {
    // ----- Ac-225 α 源：区室抽样位置 + 各向同性 + 抽样能量 -----
    G4double ekin = SampleAc225AlphaEnergy();
    G4ThreeVector pos = SampleSourcePosition();
    G4ThreeVector dir = SampleIsotropicDirection();

    fParticleGun->SetParticleDefinition(fAlpha);
    fParticleGun->SetParticleEnergy(ekin);
    fParticleGun->SetParticlePosition(pos);
    fParticleGun->SetParticleMomentumDirection(dir);

    // 诊断：打印前若干事件的抽样参数(验证区室位置 + 各向同性 + 能量)
    if (anEvent->GetEventID() < 10) {
      G4cout << "[Ac225] event " << anEvent->GetEventID()
             << "  comp=" << fCompartment
             << "  Ekin=" << ekin / MeV << " MeV"
             << "  |pos|=" << pos.mag() / um << " um"
             << "  pos(um)=" << pos / um
             << "  dir=" << dir << G4endl;
    }
  }
  else if (fSourceType == "alpha") {
    // ----- 单能 α 验证模式(任务3.1)：能量由 /gun/energy 给，位置按区室，各向同性 -----
    G4ThreeVector pos = SampleSourcePosition();
    G4ThreeVector dir = SampleIsotropicDirection();

    fParticleGun->SetParticleDefinition(fAlpha);
    fParticleGun->SetParticlePosition(pos);
    fParticleGun->SetParticleMomentumDirection(dir);

    if (anEvent->GetEventID() < 5) {
      G4cout << "[alpha] event " << anEvent->GetEventID()
             << "  comp=" << fCompartment
             << "  |pos|=" << pos.mag() / um << " um"
             << "  dir=" << dir << G4endl;
    }
  }
  else {
    // ----- proton 模式：保持原有基线行为 -----
    fParticleGun->SetParticlePosition(G4ThreeVector(fX0, fY0, fZ0));
    fParticleGun->SetParticleMomentumDirection(
      G4ThreeVector(fMomentumX, fMomentumY, fMomentumZ));
  }

  // this function is called at the begining of event
  fParticleGun->GeneratePrimaryVertex(anEvent);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double PrimaryGeneratorAction::SampleAc225AlphaEnergy() const
{
  // Ac-225 衰变链 4 个 α（能量据 NNDC，实现时已核对）：
  //   1) 225Ac -> 221Fr : 5.83 MeV  (100%)
  //   2) 221Fr -> 217At : 6.36 MeV  (100%)
  //   3) 217At -> 213Bi : 7.07 MeV  (100%)
  //   4) 经 213Bi: 97.84% -> 213Po -> 209Pb (8.40 MeV)
  //                 2.16% -> 209Tl 直接 alpha (5.87 MeV)
  // 一次衰变产生 4 个 α(各槽位各 1 个)，故 4 槽等概率抽样即可重建谱形。
  const G4double eAc = 5.83 * MeV;
  const G4double eFr = 6.36 * MeV;
  const G4double eAt = 7.07 * MeV;
  const G4double ePo = 8.40 * MeV;
  const G4double eBi = 5.87 * MeV;

  G4int slot = static_cast<G4int>(G4UniformRand() * 4.0);  // 0..3 等概率
  if (slot == 0) return eAc;
  if (slot == 1) return eFr;
  if (slot == 2) return eAt;
  return (G4UniformRand() < 0.9784) ? ePo : eBi;  // slot==3, 含 Bi-213 分支
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4ThreeVector PrimaryGeneratorAction::SampleIsotropicDirection() const
{
  // 各向同性单位向量：cosθ 在 [-1,1] 均匀，φ 在 [0,2π) 均匀
  G4double cosTheta = 2. * G4UniformRand() - 1.;
  G4double sinTheta = std::sqrt((1. - cosTheta) * (1. + cosTheta));
  G4double phi = CLHEP::twopi * G4UniformRand();
  return G4ThreeVector(sinTheta * std::cos(phi), sinTheta * std::sin(phi),
                       cosTheta);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4ThreeVector
PrimaryGeneratorAction::SampleSourcePosition() const
{
  // 按 fCompartment 在对应区室内均匀抽样源点位置
  const auto* det = static_cast<const DetectorConstruction*>(
    G4RunManager::GetRunManager()->GetUserDetectorConstruction());
  G4double Rn = det ? det->GetNucleusRadius() : 8. * um;
  G4double Rc = det ? det->GetCellRadius() : 10. * um;
  G4double Rworld = Rc + (det ? det->GetMaxRange() : 49. * um);

  if (fCompartment == "Nucleus")       return SampleInSphere(Rn);
  if (fCompartment == "Cytoplasm")     return SampleInShell(Rn, Rc);
  if (fCompartment == "Extracellular") return SampleInBoxMinusSphere(Rc, Rworld);
  // 默认 / Membrane：细胞膜面
  return SampleOnSphere(Rc);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4ThreeVector PrimaryGeneratorAction::SampleInSphere(G4double R) const
{
  // 球内均匀：r ~ R * u^(1/3)，方向各向同性
  G4double r = R * std::pow(G4UniformRand(), 1. / 3.);
  return r * SampleIsotropicDirection();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4ThreeVector
PrimaryGeneratorAction::SampleInShell(G4double Rin, G4double Rout) const
{
  // 球壳 [Rin,Rout] 内均匀：r^3 在 [Rin^3, Rout^3] 均匀
  G4double rin3 = Rin * Rin * Rin;
  G4double rout3 = Rout * Rout * Rout;
  G4double r = std::pow(rin3 + G4UniformRand() * (rout3 - rin3), 1. / 3.);
  return r * SampleIsotropicDirection();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4ThreeVector PrimaryGeneratorAction::SampleOnSphere(G4double R) const
{
  // 球面均匀 = R * 各向同性单位方向
  return R * SampleIsotropicDirection();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4ThreeVector
PrimaryGeneratorAction::SampleInBoxMinusSphere(G4double Rc, G4double wh) const
{
  // 半边长 wh 的盒内、排除半径 Rc 的球(胞外区)：拒绝采样
  for (G4int i = 0; i < 1000; ++i) {
    G4double x = (2. * G4UniformRand() - 1.) * wh;
    G4double y = (2. * G4UniformRand() - 1.) * wh;
    G4double z = (2. * G4UniformRand() - 1.) * wh;
    G4ThreeVector p(x, y, z);
    if (p.mag() > Rc) return p;
  }
  return G4ThreeVector(wh, 0., 0.);  // fallback(几乎不会到达)
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
