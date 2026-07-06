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
/// \file DetectorConstruction.cc
/// \brief Implementation of the DetectorConstruction class
///
/// 细胞几何：同心球结构 World(水盒) -> Cell(细胞膜边界) -> Nucleus(细胞核)
/// 敏感探测器(SD)置于细胞核，域采样与核打分均在核内进行。
/// 当前所有体积均为液态水（核成分后续可细化）。

#include "DetectorConstruction.hh"

#include "G4Box.hh"
#include "G4Orb.hh"
#include "G4Exception.hh"
#include "G4Material.hh"
#include "G4NistManager.hh"
#include "G4PVPlacement.hh"
#include "G4Region.hh"
#include "G4RunManager.hh"
#include "G4SDManager.hh"
#include "G4StateManager.hh"
#include "G4SystemOfUnits.hh"
#include "globals.hh"

#include "DetectorMessenger.hh"
#include "TrackerSD.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::DetectorConstruction() : G4VUserDetectorConstruction()
{
  // Create commands for interactive definition of the detector
  fDetectorMessenger = std::make_unique<DetectorMessenger>(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::~DetectorConstruction() = default;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::Construct()
{
  // Define the materials
  DefineMaterials();

  // Define volumes: 世界内含细胞，细胞内含细胞核
  G4VPhysicalVolume* World = DefineWorld();
  return World;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::DefineMaterials()
{
  // Taking water from the NIST database
  G4NistManager* nist = G4NistManager::Instance();
  if (!fMat) {
    fMat = nist->FindOrBuildMaterial("G4_WATER");
  }

  // Print the material information
  G4cout << "Material: " << fMat->GetName() << G4endl;
  G4cout << *(G4Material::GetMaterialTable()) << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::CheckConsistency()
{
  if (!(fMaxRange > 0.)) {
    G4ExceptionDescription msg;
    msg << "fMaxRange must be > 0.\n"
        << "Please consider using the formula by Tabata for the calculation "
        << "of the maximum range of secondary electrons:\n"
        << "\t https://doi.org/10.1016/0029-554X(72)90463-6" << G4endl
        << "Note: 此参数决定世界体积在细胞外的水层厚度（次级电子平衡）。"
        << G4endl;
    G4Exception("DetectorConstruction::CheckConsistency()", "DetCons0001",
                FatalException, msg);
  }
  if (!(fNucleusRadius > 0.)) {
    G4ExceptionDescription msg;
    msg << "fNucleusRadius must be > 0.\n"
        << "细胞核半径必须为正，否则无法在核内放置采样位点。" << G4endl;
    G4Exception("DetectorConstruction::CheckConsistency()", "DetCons0002",
                FatalException, msg);
  }
  if (!(fCellRadius > fNucleusRadius)) {
    G4ExceptionDescription msg;
    msg << "fCellRadius must be > fNucleusRadius.\n"
        << "细胞半径必须大于细胞核半径（同心球结构）。" << G4endl;
    G4Exception("DetectorConstruction::CheckConsistency()", "DetCons0003",
                FatalException, msg);
  }
  if (!(fSiteRadius > 0. && fSiteRadius <= fNucleusRadius)) {
    G4ExceptionDescription msg;
    msg << "fSiteRadius must be > 0 and <= fNucleusRadius.\n"
        << "域(site)半径必须为正且不大于细胞核半径。" << G4endl;
    G4Exception("DetectorConstruction::CheckConsistency()", "DetCons0004",
                FatalException, msg);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::DefineWorld()
{
  // Ensure geometry consistency
  CheckConsistency();

  // 世界为水盒，半边长 = 细胞半径 + 次级电子最大射程
  // 保证细胞外有足够水层维持次级电子平衡
  G4double WorldHalf = fCellRadius + fMaxRange;

  auto* sWorld = new G4Box("World", WorldHalf, WorldHalf, WorldHalf);
  auto* lWorld = new G4LogicalVolume(sWorld, fMat, "World", 0, 0, 0);
  G4VPhysicalVolume* World =
    new G4PVPlacement(0, G4ThreeVector(), lWorld, "World", nullptr, false, 0);

  // 细胞球置于世界内（位于原点）
  DefineCell(World);

  // Print the world volume information
  PrintParameters(World);

  return World;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume*
DetectorConstruction::DefineCell(G4VPhysicalVolume* mother)
{
  // 细胞(膜边界)球，材料为水（胞质近似），中心位于原点
  auto* sCell = new G4Orb("Cell", fCellRadius);
  auto* lCell = new G4LogicalVolume(sCell, fMat, "Cell", 0, 0, 0);
  G4VPhysicalVolume* Cell = new G4PVPlacement(
    0, G4ThreeVector(), lCell, "Cell", mother->GetLogicalVolume(), false, 0);

  // 细胞核球置于细胞内（位于原点）
  DefineNucleus(Cell);

  // Print the cell volume information
  PrintParameters(Cell);

  return Cell;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume*
DetectorConstruction::DefineNucleus(G4VPhysicalVolume* mother)
{
  // 细胞核球，材料为水（核成分暂以水近似），中心位于原点
  auto* sNucleus = new G4Orb("Nucleus", fNucleusRadius);
  auto* lNucleus = new G4LogicalVolume(sNucleus, fMat, "Nucleus", 0, 0, 0);
  G4VPhysicalVolume* Nucleus =
    new G4PVPlacement(0, G4ThreeVector(), lNucleus, "Nucleus",
                      mother->GetLogicalVolume(), false, 0);

  // Print the nucleus volume information
  PrintParameters(Nucleus);

  return Nucleus;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::ConstructSDandField()
{
  // 敏感探测器置于细胞核：域级采样与核级打分均在核内进行
  G4String SDname = "Nucleus";

  auto* aSD = new TrackerSD(SDname, "TrackerHitColl");

  // Register the SD with the SD manager
  G4SDManager* SDMan = G4SDManager::GetSDMpointer();
  SDMan->AddNewDetector(aSD);

  // Set the SD to the Nucleus logical volume
  SetSensitiveDetector("Nucleus", aSD, true);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetMaterial(const G4String& name)
{
  auto* nist = G4NistManager::Instance();

  G4Material* mat = nist->FindOrBuildMaterial(name);

  if (!mat) {
    G4ExceptionDescription ed;
    ed << "Material '" << name << "' not found. Falling back to G4_WATER.";
    G4Exception("DetectorConstruction::SetMaterial", "MTK001", JustWarning, ed);
    mat = nist->FindOrBuildMaterial("G4_WATER");
  }

  if (fMat != mat) {
    fMat = mat;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::PrintParameters(G4VPhysicalVolume* physVol) const
{
  if (!physVol) {
    G4cerr << "Error: Physical volume is null!" << G4endl;
    return;
  }

  G4ThreeVector pos = physVol->GetObjectTranslation();
  G4LogicalVolume* logVol = physVol->GetLogicalVolume();
  G4VSolid* solid = logVol->GetSolid();
  G4Material* material = logVol->GetMaterial();

  G4cout << "\n================ Volume Parameters ================" << G4endl;
  G4cout << "Physical Volume Name: " << physVol->GetName() << G4endl;
  G4cout << "Position: " << pos / mm << " mm" << G4endl;
  G4cout << "Material: " << material->GetName() << G4endl;

  // Detect shape and print dimensions accordingly
  if (auto box = dynamic_cast<G4Box*>(solid)) {
    G4cout << "Shape: Box" << G4endl;
    G4cout << "Dimensions (full): " << 2 * box->GetXHalfLength() / mm << " x "
           << 2 * box->GetYHalfLength() / mm << " x "
           << 2 * box->GetZHalfLength() / mm << " mm" << G4endl;
  }
  else if (auto orb = dynamic_cast<G4Orb*>(solid)) {
    G4cout << "Shape: Orb (sphere)" << G4endl;
    G4cout << "Radius: " << orb->GetRadius() / um << " um" << G4endl;
  }
  else {
    G4cout << "Shape: Unknown or not handled yet." << G4endl;
  }

  G4cout << "===================================================\n" << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
