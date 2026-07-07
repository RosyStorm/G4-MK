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
/// \file PrimaryGeneratorAction.hh
/// \brief Definition of the PrimaryGeneratorAction class

#ifndef PrimaryGeneratorAction_h
#define PrimaryGeneratorAction_h 1

#include "G4ParticleTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4ThreeVector.hh"
#include "G4VUserPrimaryGeneratorAction.hh"
#include "globals.hh"

#include <memory>

class PrimaryGeneratorMessenger;
class G4ParticleGun;
class G4ParticleDefinition;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
  public:
    PrimaryGeneratorAction();
    virtual ~PrimaryGeneratorAction() override;
    void GeneratePrimaries(G4Event*) override;

    G4ParticleGun* GetParticleGun() { return fParticleGun.get(); }

    // Setters for the particle source
    void SetPositionZ(G4double z) { fZ0 = z; }
    void SetSourceType(const G4String& t) { fSourceType = t; }
    G4String GetSourceType() const { return fSourceType; }
    void SetCompartment(const G4String& c) { fCompartment = c; }
    G4String GetCompartment() const { return fCompartment; }

  private:
    // Ac-225 α 源辅助方法
    G4double SampleAc225AlphaEnergy() const;            // 抽样 Ac-225 衰变链 α 动能
    G4ThreeVector SampleIsotropicDirection() const;     // 各向同性单位方向
    G4ThreeVector SampleSourcePosition() const;         // 按 fCompartment 抽样源点位置
    G4ThreeVector SampleInSphere(G4double R) const;                  // 球内均匀
    G4ThreeVector SampleInShell(G4double Rin, G4double Rout) const;  // 球壳内均匀
    G4ThreeVector SampleOnSphere(G4double R) const;                  // 球面均匀
    G4ThreeVector SampleInBoxMinusSphere(G4double Rc, G4double wh) const; // 盒内排除球(胞外)

    std::unique_ptr<G4ParticleGun> fParticleGun;
    std::unique_ptr<PrimaryGeneratorMessenger> fGunMessenger;

    // Source properties
    G4ParticleDefinition* fParticle = nullptr;   // 质子(基线对比用)
    G4ParticleDefinition* fAlpha = nullptr;      // α 粒子
    G4String fSourceType = "ac225";              // proton | ac225
    G4String fCompartment = "Membrane";          // Nucleus | Cytoplasm | Membrane | Extracellular
    G4double fX0 = 0.;
    G4double fY0 = 0.;
    G4double fZ0 = -10.2 * um;
    G4double fMomentumX = 0.;
    G4double fMomentumY = 0.;
    G4double fMomentumZ = 1.;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
