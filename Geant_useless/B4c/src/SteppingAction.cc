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
/// \file electromagnetic/TestEm3/src/SteppingAction.cc
/// \brief Implementation of the SteppingAction class
//
// $Id: SteppingAction.cc 98762 2016-08-09 14:08:07Z gcosmo $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "SteppingAction.hh"

#include "B4cDetectorConstruction.hh"
#include "B4cEventAction.hh"

#include "G4Positron.hh"
#include "G4RunManager.hh"
#include "G4PhysicalConstants.hh"
#include "read_simu_tree.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SteppingAction::SteppingAction(B4cEventAction* evt)
:G4UserSteppingAction(),fEventAct(evt) 
{
  fCurrTrack = 0;
  fDsY = 0;
  fDsZ = 0;
  fDsYZ = 0;
  fDsZY = 0;
  fDsXY = 0;
  fDsXZ = 0;
  fdEY = 0;
  fdEZ = 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SteppingAction::~SteppingAction()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SteppingAction::UserSteppingAction(const G4Step* aStep)
{
  //track informations
  const G4StepPoint* prePoint = aStep->GetPreStepPoint();   
  const G4StepPoint* endPoint = aStep->GetPostStepPoint();
  const G4ParticleDefinition* particle = aStep->GetTrack()->GetDefinition(); 
    
  //if World, return
  //
  G4VPhysicalVolume* volume = prePoint->GetTouchableHandle()->GetVolume();    
  //if sum of absorbers do not fill exactly a layer: check material, not volume.
  if (volume == G4LogicalVolumeStore::GetInstance()->GetVolume("world");) return;
 
  //here we are in an absorber. Locate it
  G4int absorNum  = prePoint->GetTouchableHandle()->GetCopyNumber(0);
  
  // collect energy deposit taking into account track weight
  G4double edep = aStep->GetTotalEnergyDeposit()*aStep->GetTrack()->GetWeight();
  
  // collect step length of charged particles
//  G4double stepl = 0.;
//  if (particle->GetPDGCharge() != 0.) {
//    stepl = aStep->GetStepLength();
//  }
  
  //  G4cout << "Nabs= " << absorNum << "   edep(keV)= " << edep << G4endl;
  
  
  //longitudinal profile of edep per absorber
  if (edep>0.) {
    G4AnalysisManager::Instance()->FillH1(kMaxAbsor+absorNum, 
                                          G4double(layerNum+1), edep);
  }

////  example of Birk attenuation
  G4double destep   = aStep->GetTotalEnergyDeposit();
  G4double response = BirksAttenuation(aStep);
  G4ThreeVector prepos = prePoint->GetPosition();
  G4ThreeVector postpos = endPoint->GetPosition();
  
  if(fDsY > 0.3125){
    myhit_view.push_back(0);
    myhit_y.push_back((postpos[0]-prepos[0])*cm/2.);
    myhit_z.push_back((postpos[1]-prepos[1])*cm/2.);
    myhit_x.push_back((postpos[2]-prepos[2])*cm/2.);
    myhit_dE.push_back(fdEY);
    myhit_ds.push_back(TMath::Sqrt(fDsY**2 + fDsYZ**2 + fDsXY**2));
    fDsY = 0;
    fDsYZ = 0;
    fDsXY = 0;
    fdEY = 0;
  }
  if(fDsZ > 0.3125){
    myhit_view.push_back(1);
    myhit_y.push_back((postpos[0]-prepos[0])*cm/2.);
    myhit_z.push_back((postpos[1]-prepos[1])*cm/2.);
    myhit_x.push_back((postpos[2]-prepos[2])*cm/2.);
    myhit_dE.push_back(fdEZ);
    myhit_ds.push_back(TMath::Sqrt(fDsY**2 + fDsZY**2 + fDsXZ**2));
    fDsZ = 0;
    fDsZY = 0;
    fDsXZ = 0;
    fdEZ = 0;
  }
  if(!fCurrTrack){
    fCurrTrack = aStep->GetTrack();
    preposY = ((int)prepos[0]*cm/0.3125) + ( prepos[0]*cm - ((int)prepos[0]*cm/0.3125) );
    preposZ = ((int)prepos[1]*cm/0.3125) + ( prepos[1]*cm - ((int)prepos[1]*cm/0.3125) );
    fDsY = postpos[0]*cm - preposY;
    fDsZ = postpos[1]*cm - preposZ;
    fDsXY = postpos[2]*cm - postpos[2]*cm;
    fDsXZ = postpos[2]*cm - postpos[2]*cm;
    fdEY = response/2.;
    fdEZ = response/2.;
  }
  else if (fCurrTrack == aStep->GetTrack()){
    fDsY += postpos[0]*cm - prepos[0]*cm;
    fDsZ += postpos[1]*cm - prepos[1]*cm;
    fDsXY += postpos[2]*cm - prepos[2]*cm;
    fDsXZ += postpos[2]*cm - prepos[2]*cm;
    fdEY += resonse/2.;
    fdEZ += resonse/2.;
  }
  else{
    fCurrTrack = aStep->GetTrack();
    preposY = ((int)prepos[0]*cm/0.3125) + ( prepos[0]*cm - ((int)prepos[0]*cm/0.3125) );
    preposZ = ((int)prepos[1]*cm/0.3125) + ( prepos[1]*cm - ((int)prepos[1]*cm/0.3125) );
    fDsY = postpos[0]*cm - preposY;
    fDsZ = postpos[1]*cm - preposZ;
    fDsXY = postpos[2]*cm - postpos[2]*cm;
    fDsXZ = postpos[2]*cm - postpos[2]*cm;
    fdEY = response/2.;
    fdEZ = response/2.;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double SteppingAction::BirksAttenuation(const G4Step* aStep)
{
  //Example of Birk attenuation law in organic scintillators.
  //adapted from Geant3 PHYS337. See MIN 80 (1970) 239-244
  //
  G4Material* material = aStep->GetTrack()->GetMaterial();
  G4double birk1       = material->GetIonisation()->GetBirksConstant();
  G4double destep      = aStep->GetTotalEnergyDeposit();
  G4double stepl       = aStep->GetStepLength();  
  G4double charge      = aStep->GetTrack()->GetDefinition()->GetPDGCharge();
  cout << "Birk's constant : " << birk1 << endl;
  //
  G4double response = destep;
  if (birk1*destep*stepl*charge != 0.)
  {
    response = destep/(1. + birk1*destep/stepl);
  }
 return response;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

