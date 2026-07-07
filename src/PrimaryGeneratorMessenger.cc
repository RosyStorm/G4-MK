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
/// \file PrimaryGeneratorMessenger.cc
/// \brief Implementation of the PrimaryGeneratorMessenger class

#include "PrimaryGeneratorMessenger.hh"

#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIdirectory.hh"

#include "PrimaryGeneratorAction.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorMessenger::PrimaryGeneratorMessenger(
  PrimaryGeneratorAction* primaryGen)
  : G4UImessenger(), fPrimaryAction(primaryGen)
{
  fGunDir = std::make_unique<G4UIdirectory>("/beam/");
  fGunDir->SetGuidance("PrimaryGenerator control");

  fZ0Cmd =
    std::make_unique<G4UIcmdWithADoubleAndUnit>("/beam/position/Z0", this);
  fZ0Cmd->SetGuidance("Set Z coordinate of the particle source");
  fZ0Cmd->SetParameterName("posZ", false);
  fZ0Cmd->SetDefaultValue(0.0);
  fZ0Cmd->SetDefaultUnit("cm");
  fZ0Cmd->SetUnitCategory("Length");
  fZ0Cmd->AvailableForStates(G4State_PreInit, G4State_Idle);

  fSourceDir = std::make_unique<G4UIdirectory>("/source/");
  fSourceDir->SetGuidance("Particle source type control");

  fSourceTypeCmd = std::make_unique<G4UIcmdWithAString>("/source/type", this);
  fSourceTypeCmd->SetGuidance(
    "Source type: proton (baseline) | ac225 (Ac-225 alpha, compartment set by /source/compartment)"
    " | alpha (mono-energetic alpha for range validation, energy from /gun/energy)");
  fSourceTypeCmd->SetParameterName("type", false);
  fSourceTypeCmd->SetCandidates("proton ac225 alpha");
  fSourceTypeCmd->SetDefaultValue("ac225");
  fSourceTypeCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

  fCompartmentCmd = std::make_unique<G4UIcmdWithAString>("/source/compartment", this);
  fCompartmentCmd->SetGuidance(
    "Source compartment (ac225): Nucleus | Cytoplasm | Membrane | Extracellular");
  fCompartmentCmd->SetParameterName("comp", false);
  fCompartmentCmd->SetCandidates("Nucleus Cytoplasm Membrane Extracellular");
  fCompartmentCmd->SetDefaultValue("Membrane");
  fCompartmentCmd->AvailableForStates(G4State_PreInit, G4State_Idle);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorMessenger::~PrimaryGeneratorMessenger() = default;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PrimaryGeneratorMessenger::SetNewValue(G4UIcommand* command,
                                            G4String newValue)
{
  if (command == fZ0Cmd.get()) {
    fPrimaryAction->SetPositionZ(fZ0Cmd->GetNewDoubleValue(newValue));
  }
  if (command == fSourceTypeCmd.get()) {
    fPrimaryAction->SetSourceType(newValue);
  }
  if (command == fCompartmentCmd.get()) {
    fPrimaryAction->SetCompartment(newValue);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......