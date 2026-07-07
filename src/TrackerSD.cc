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
/// \file TrackerSD.cc
/// \brief Implementation of the TrackerSD class

#include "TrackerSD.hh"

#include "G4AnalysisManager.hh"
#include "G4Event.hh"
#include "G4Exception.hh"
#include "G4ParticleDefinition.hh"
#include "G4PrimaryParticle.hh"
#include "G4PrimaryVertex.hh"
#include "G4ProcessType.hh"
#include "G4RunManager.hh"
#include "G4SDManager.hh"
#include "G4Step.hh"
#include "G4StepPoint.hh"
#include "G4SystemOfUnits.hh"
#include "G4VProcess.hh"
#include "Randomize.hh"
#include "globals.hh"

#include "DetectorConstruction.hh"
#include "PrimaryGeneratorAction.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
TrackerSD::TrackerSD(const G4String& sdName, const G4String& hitsCollectionName,
                     const int /*depthIndex*/)
  : G4VSensitiveDetector(sdName)
{
  collectionName.insert(hitsCollectionName);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

TrackerSD::~TrackerSD() = default;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void TrackerSD::Initialize(G4HCofThisEvent* hce)
{
  // Create hits collection
  fHitsCollection =
    new TrackerHitColl(SensitiveDetectorName, collectionName[0]);

  // Add this collection in hce
  G4int collID =
    G4SDManager::GetSDMpointer()->GetCollectionID(collectionName[0]);

  hce->AddHitsCollection(collID, fHitsCollection);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool TrackerSD::ProcessHits(G4Step* aStep, G4TouchableHistory*)
{
  // Get the track
  G4Track* aTrack = aStep->GetTrack();

  // Get Pre and Post step points
  G4StepPoint* preStepPoint = aStep->GetPreStepPoint();
  G4StepPoint* postStepPoint = aStep->GetPostStepPoint();

  // Get ParentID
  G4int parentID = aTrack->GetParentID();

  // Only for primary particles: 记录核边界入射/出射动能
  // 注: 用 fGeomBoundary 状态位判定。对 DNA 物理下的 α 离子边界步该状态位可能不触发,
  //     故 α 的 T_in/T_out 可能为空(诊断量, 不影响 z/y/ε 打分)。
  //     若任务6.2 需精确能量平衡, 改用 G4UserSteppingAction 按体积过渡判定。
  if (parentID == 0) {
    if (preStepPoint->GetStepStatus() == fGeomBoundary) {
      G4AnalysisManager::Instance()->FillH1(12, preStepPoint->GetKineticEnergy());
    }
    if (postStepPoint->GetStepStatus() == fGeomBoundary) {
      G4AnalysisManager::Instance()->FillH1(13, postStepPoint->GetKineticEnergy());
    }
  }

  // Get the energy deposit
  G4double edep = aStep->GetTotalEnergyDeposit();
  if (edep == 0.) return false;

  // Declare a new hit
  auto newHit = new TrackerHit();

  newHit->SetEdep(edep);

  // Get the position of the hit
  G4ThreeVector Pos = postStepPoint->GetPosition();
  newHit->SetPosition(Pos);

  // Add the hit to the collection
  fHitsCollection->insert(newHit);

  return true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void TrackerSD::EndOfEvent(G4HCofThisEvent*)
{
  const auto nofHits = fHitsCollection->entries();

  if (verboseLevel > 0) {
    G4cout << "  Hits Collection: in this event, we have " << nofHits
           << " hits in the tracker." << G4endl;
  }

  // 注：miss 事件(nofHits==0, 如膜/胞外源 α 未命中核)不再早返回——
  // ntuple 用 hitFlag=0 记录(任务4.2)；直方图 f_{n,1} 由各 if(nucleusEdep>0)/if(nHint>0) 自然过滤。

  // Get the DetectorConstruction instance
  const DetectorConstruction* detConstruction =
    static_cast<const DetectorConstruction*>(
      G4RunManager::GetRunManager()->GetUserDetectorConstruction());

  // Declare variables for the dimensions of the detector
  G4double nucleusRadius = 0.;
  G4double site_radius = 0.;
  G4double DetDensity = 0.;
  G4ThreeVector hitPos;
  G4ThreeVector randCenterPos;

  if (detConstruction) {
    site_radius = detConstruction->GetSiteRadius();
    nucleusRadius = detConstruction->GetNucleusRadius();
    auto mat = G4Material::GetMaterial(detConstruction->GetMaterial());
    if (mat)
      DetDensity = mat->GetDensity();
    else {
      G4ExceptionDescription msg;
      msg << "Material not found."
          << "Something unexpected has occurred." << G4endl;
      G4Exception("TrackerSD::EndOfEvent()", "TrackerSD002", JustWarning, msg);
    }
  }
  else {
    G4ExceptionDescription msg;
    msg << "Detector construction not found."
        << "Something unexpected has occurred." << G4endl;
    G4Exception("TrackerSD::EndOfEvent()", "TrackerSD001", JustWarning, msg);
  }
  G4int nHsel = 0;  // number of hits in the selection region
  G4int nHsite = 0;  // number of hits in the site
  G4int nHint = 0;  // number of hits both in site and selection region
  G4double evtEdep = 0.;  // 域内能量沉积(单事件)
  G4double nucleusEdep = 0.;  // 核内总能量沉积(单事件, 任务4.1 用于 z_n)

  if (nofHits > 0) {
    const G4int maxTries = 1000;
    G4int tries = 0;
    G4bool found = false;

    while (tries < maxTries && !found) {
      std::size_t randHit = static_cast<std::size_t>(G4UniformRand() * nofHits);
      auto hit = (*fHitsCollection)[randHit];
      G4ThreeVector hitPosition = hit->GetPosition();

      // 判定该 hit 是否落在细胞核内（核球为 hit 选择区）
      if (hitPosition.mag() < nucleusRadius)
      {
        hitPos = hitPosition;  // store valid hit position

        found = true;
      }
      ++tries;
    }

    if (found) {
      G4double site_radius2 = site_radius * site_radius;

      // Random placement of a sphere center containing the hit
      G4double xRand, yRand, zRand, randRad2;
      do {
        xRand = (2 * G4UniformRand() - 1) * site_radius;
        yRand = (2 * G4UniformRand() - 1) * site_radius;
        zRand = (2 * G4UniformRand() - 1) * site_radius;
        randRad2 = xRand * xRand + yRand * yRand + zRand * zRand;
      } while (randRad2 > site_radius2);

      randCenterPos = G4ThreeVector(xRand + hitPos.x(), yRand + hitPos.y(),
                                    zRand + hitPos.z());
    }

    else {
      G4cout
        << "In this event, no hits were found in the nucleus "
        << "after trying " << tries << " times." << G4endl
        << "Please consider checking the source position / nucleus radius."
        << G4endl;
      return;  // Skip analysis if no valid hit was found
    }
  }

  // Loop over hits
  for (std::size_t jj = 0; jj < nofHits; jj++) {
    auto hit2 = (*fHitsCollection)[jj];
    G4ThreeVector hit2Position = hit2->GetPosition();

    // 核总能量累加(所有 hit 均在核内, SD 在核上)
    nucleusEdep += hit2->GetEdep();

    // 判定该 hit 是否落在细胞核内
    G4bool inNucleus = (hit2Position.mag() < nucleusRadius);

    // Check if the hit is within the site
    G4double dist = (hit2Position - randCenterPos).mag();
    if (dist < site_radius) {
      nHsite++;
      evtEdep += hit2->GetEdep();
    }

    // Check if the hit is within the nucleus
    if (inNucleus) nHsel++;

    // Check if the hit is within the nucleus and the site
    if (dist < site_radius && inNucleus) nHint++;
  }

  // Access the analysis manager
  auto analysisManager = G4AnalysisManager::Instance();

  // ===== 核级比能 z_n (任务4.1)：核内总沉积 / 核质量，无加权(核是整个位点) =====
  G4double massNucleus =
    (4. / 3.) * CLHEP::pi * nucleusRadius * nucleusRadius * nucleusRadius
    * DetDensity;
  G4double z_n = nucleusEdep / massNucleus;
  if (nucleusEdep > 0.) {
    G4double zn = z_n / gray;
    G4int idF = analysisManager->GetH1Id("fzn");
    G4int idZF = analysisManager->GetH1Id("znfzn");
    G4int idZ2F = analysisManager->GetH1Id("z2nfzn");
    if (idF >= 0) analysisManager->FillH1(idF, zn, 1.);            // → z̄_{n,F}
    if (idZF >= 0) analysisManager->FillH1(idZF, zn, zn);          // → z̄_{n,D}
    if (idZ2F >= 0) analysisManager->FillH1(idZ2F, zn, zn * zn);   // 高阶矩
  }

  // Calculate lineal energy
  G4double y = (evtEdep) / ((4. / 3.) * site_radius);  // in keV/um

  // Calculate specific energy
  G4double mass = ((4. / 3.) * CLHEP::pi * site_radius * site_radius
                   * site_radius * DetDensity);
  G4double z = (evtEdep / mass);

  // Fill histograms
  if (nHint > 0) {
    // Define the weight
    G4double wght = G4double(nHsel) / G4double(nHint);

    // Histogram 0: Single-event energy imparted in keV
    analysisManager->FillH1(0, evtEdep / keV, wght);

    // Histogram 1: Weighted single-event energy imparted in keV
    analysisManager->FillH1(1, evtEdep / keV, (evtEdep * wght / keV));

    // Histogram 2: Squared weighted energy imparted per event in keV^2
    analysisManager->FillH1(2, evtEdep / keV,
                            ((evtEdep * evtEdep * wght) / (keV * keV)));

    // Histogram 3: Lineal energy in keV/um
    analysisManager->FillH1(3, y / (keV / um), wght);

    // Histogram 4: Dose-weighted lineal energy in keV/um
    analysisManager->FillH1(4, y / (keV / um), (y * wght / (keV / um)));

    // Histogram 5: Squared weighted lineal energy in (keV/um)^2
    analysisManager->FillH1(5, y / (keV / um),
                            (y * y * wght / ((keV / um) * (keV / um))));

    // Histogram 6: Specific energy in keV/g
    analysisManager->FillH1(6, z / gray, wght);

    // Histogram 7: Weighted specific energy in keV/g
    analysisManager->FillH1(7, z / gray, (z * wght / gray));

    // Histogram 8: Squared-weighted specific energy in (keV/g)^2
    analysisManager->FillH1(8, z / gray, (z * z * wght / (gray * gray)));

    // Histogram 9: Number of hits in the selection region
    analysisManager->FillH1(9, nHsel, 1.);

    // Histogram 10: Number of hits in the site
    analysisManager->FillH1(10, nHsite, 1.);

    // Histogram 11: Number of hits both in site and selection region
    analysisManager->FillH1(11, nHint, 1.);

    // Histogram 14: Number of hits in site vs Energy imparted per event
    analysisManager->FillH2(0, evtEdep / keV, nHsite, 1.);
  }
  // (else: miss 事件 nofHits==0 —— 对膜/胞外源属正常, 不再警告; ntuple 用 hitFlag=0 记录)

  // ===== 任务4.2：逐事件配对 ntuple (hit & miss 都填一行) =====
  // 取初级 α 动能(从 primary vertex)与区室编号(从 PrimaryGeneratorAction)
  G4double alphaE = 0.;
  const G4Event* evt = G4RunManager::GetRunManager()->GetCurrentEvent();
  if (evt && evt->GetNumberOfPrimaryVertex() > 0) {
    const G4PrimaryVertex* pv = evt->GetPrimaryVertex(0);
    if (pv) {
      const G4PrimaryParticle* pp = pv->GetPrimary(0);
      if (pp) alphaE = pp->GetKineticEnergy();
    }
  }
  G4int compId = 4;  // 未知
  const auto* pga = dynamic_cast<const PrimaryGeneratorAction*>(
    G4RunManager::GetRunManager()->GetUserPrimaryGeneratorAction());
  if (pga) {
    const G4String c = pga->GetCompartment();
    if (c == "Nucleus") compId = 0;
    else if (c == "Cytoplasm") compId = 1;
    else if (c == "Membrane") compId = 2;
    else if (c == "Extracellular") compId = 3;
  }
  G4double wght = (nHint > 0) ? G4double(nHsel) / G4double(nHint) : 0.;
  analysisManager->FillNtupleIColumn(0, evt ? evt->GetEventID() : 0);
  analysisManager->FillNtupleDColumn(1, alphaE / MeV);
  analysisManager->FillNtupleDColumn(2, evtEdep / keV);
  analysisManager->FillNtupleDColumn(3, z / gray);          // z_d
  analysisManager->FillNtupleDColumn(4, nucleusEdep / keV);
  analysisManager->FillNtupleDColumn(5, z_n / gray);        // z_n
  analysisManager->FillNtupleDColumn(6, wght);
  analysisManager->FillNtupleIColumn(7, nHsel);
  analysisManager->FillNtupleIColumn(8, nHsite);
  analysisManager->FillNtupleIColumn(9, nHint);
  analysisManager->FillNtupleIColumn(10, nofHits > 0 ? 1 : 0);  // hitFlag
  analysisManager->FillNtupleIColumn(11, compId);               // compartment
  analysisManager->AddNtupleRow();
}