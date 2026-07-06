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
/// \file DetectorConstruction.hh
/// \brief Definition of the DetectorConstruction class

#ifndef DetectorConstruction_H
#define DetectorConstruction_H 1

#include "G4Material.hh"
#include "G4SystemOfUnits.hh"
#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"

#include <memory>

class DetectorMessenger;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class DetectorConstruction : public G4VUserDetectorConstruction
{
  public:
    // Constructor and destructor
    DetectorConstruction();
    ~DetectorConstruction() override;

    G4VPhysicalVolume* Construct() override;
    void ConstructSDandField() override;

    // SET Methods
    void SetMaterial(const G4String& name);
    void SetMaxRange(const G4double& range) { fMaxRange = range; }
    void SetCellRadius(const G4double& r) { fCellRadius = r; }
    void SetNucleusRadius(const G4double& r) { fNucleusRadius = r; }
    void SetSiteRadius(const G4double& siteRadius) { fSiteRadius = siteRadius; }
    void PrintParameters(G4VPhysicalVolume*) const;
    void CheckConsistency();

    // GET Methods
    G4String GetMaterial() const
    {
      return fMat ? fMat->GetName() : G4String("undefined");
    }
    G4double GetMaxRange() const { return fMaxRange; }
    G4double GetCellRadius() const { return fCellRadius; }
    G4double GetNucleusRadius() const { return fNucleusRadius; }
    G4double GetSiteRadius() const { return fSiteRadius; }

  private:
    void DefineMaterials();
    G4VPhysicalVolume* DefineWorld();
    G4VPhysicalVolume* DefineCell(G4VPhysicalVolume* mother);
    G4VPhysicalVolume* DefineNucleus(G4VPhysicalVolume* mother);
    std::unique_ptr<DetectorMessenger> fDetectorMessenger;

    G4Material* fMat = nullptr;
    G4double fMaxRange = 8.25 * um;
    G4double fCellRadius = 10. * um;      // 细胞(膜)半径 R_cell
    G4double fNucleusRadius = 8. * um;    // 细胞核半径 R_n
    G4double fSiteRadius = 0.5 * um;      // 域(site)半径 r_d
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
