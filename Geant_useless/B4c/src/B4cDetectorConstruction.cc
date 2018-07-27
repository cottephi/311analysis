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
// $Id: B4cDetectorConstruction.cc 101905 2016-12-07 11:34:39Z gunter $
// 
/// \file B4cDetectorConstruction.cc
/// \brief Implementation of the B4cDetectorConstruction class

#include "B4cDetectorConstruction.hh"
#include "B4cCalorimeterSD.hh"
#include "G4Material.hh"
#include "G4NistManager.hh"

#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVReplica.hh"
#include "G4GlobalMagFieldMessenger.hh"
#include "G4AutoDelete.hh"

#include "G4SDManager.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4ThreadLocal 
G4GlobalMagFieldMessenger* B4cDetectorConstruction::fMagFieldMessenger = 0; 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B4cDetectorConstruction::B4cDetectorConstruction()
 : G4VUserDetectorConstruction(),
   fCheckOverlaps(false),
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B4cDetectorConstruction::~B4cDetectorConstruction()
{ 
}  

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* B4cDetectorConstruction::Construct()
{
  // Define materials 
  DefineMaterials();
  
  // Define volumes
  return DefineVolumes();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B4cDetectorConstruction::DefineMaterials()
{ 
  // Lead material defined using NIST Manager
  auto nistManager = G4NistManager::Instance();
  
  // Liquid argon material
  G4double a;  // mass of a mole;
  G4double z;  // z=mean number of protons;  
  G4double density; 
  new G4Material("liquidArgon", z=18., a= 39.95*g/mole, density= 1.390*g/cm3);
         // The argon by NIST Manager is a gas with a different density

  // Print materials
  G4cout << *(G4Material::GetMaterialTable()) << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* B4cDetectorConstruction::DefineVolumes()
{
  // Geometry parameters
  G4double Zlength = 3.*m;
  G4double Xdrift =  1.*m;
  G4double Ywidth  = 1.*m;

  // Get materials
  auto LAr = G4Material::GetMaterial("liquidArgon");
  
  if ( ! LAr ) {
    G4ExceptionDescription msg;
    msg << "Cannot retrieve materials already defined."; 
    G4Exception("B4DetectorConstruction::DefineVolumes()",
      "MyCode0001", FatalException, msg);
  }  
   
  //     
  // World
  //
  
  auto world_S 
    = new G4Box("world",           // its name
                 Ywidth*2, Zlength*2, Xdrift*2); // its size
                         
  auto world_LV
    = new G4LogicalVolume(
                 world_S,           // its solid
                 LAr,  // its material
                 "world");         // its name
                                   
  auto world_PV
    = new G4PVPlacement(
                 0,                // no rotation
                 G4ThreeVector(0,Zlength,0),  // at (0,0,0)
                 world_LV,          // its logical volume                         
                 "world",          // its name
                 0,                // its mother  volume
                 false,            // no boolean operation
                 0,                // copy number
                 false);  // checking overlaps 
  //                                        
  // Visualization attributes
  //
  auto invisibility= new G4VisAttributes(G4Colour(0,0,0,0));
  invisibility->SetVisibility(true);
  world_LV->SetVisAttributes(invisibility);
  
  auto FiVol_S 
    = new G4Box("FiducialVolume_S",           // its name
                 Ywidth, Zlength, Xdrift); // its size
                         
  auto FiVol_LV
    = new G4LogicalVolume(
                 FiVol_S,           // its solid
                 LAr,  // its material
                 "FiducialVolume_LV");         // its name
                                   
  auto FiVol_PV
    = new G4PVPlacement(
                 0,                // no rotation
                 G4ThreeVector(0,Zlength/2,0),  // at (0,0,0)
                 FiVol_LV,          // its logical volume                         
                 "FiducialVolume_PV",          // its name
                 world_LV,                // its mother  volume
                 false,            // no boolean operation
                 0,                // copy number
                 false);  // checking overlaps 
  //                                        
  // Visualization attributes
  //
  auto simpleBoxVisAtt= new G4VisAttributes(G4Colour(0,1.0,0,.5));
  simpleBoxVisAtt->SetVisibility(true);
  FiVol_LV->SetVisAttributes(simpleBoxVisAtt);
  
  G4ElectricField* fEMfield = new G4UniformElectricField( G4ThreeVector(0., 0., -500.*volt/cm) );

  // Always return the physical World
  //
  return world_PV;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B4cDetectorConstruction::ConstructSDandField(){
  // Sensitive detectors
  auto LArSD 
    = new B4cCalorimeterSD("LArSD", "LArHitsCollection", 0);
  G4SDManager::GetSDMpointer()->AddNewDetector(LArSD);
  SetSensitiveDetector("FiducialVolume_LV",absoSD);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
