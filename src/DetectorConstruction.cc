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
/// \file electromagnetic/TestEm14/src/DetectorConstruction.cc
/// \brief Implementation of the DetectorConstruction class
//
// $Id: DetectorConstruction.cc 68313 2013-03-21 18:15:21Z maire $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//探测器构建函数
#include "DetectorConstruction.hh"
#include "DetectorMessenger.hh"

#include "G4Material.hh"
#include "G4Box.hh"
#include "G4Orb.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"

#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"

#include "G4NistManager.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::DetectorConstruction()
:G4VUserDetectorConstruction(),physiWorld(0), targetMaterial(0),
 fDetectorMessenger(0)
{
    worldR =        3.*mm;
    targetsizex = 1.0*mm;
    targetsizey = 1.0*mm;
    targetsizez = 10.0*um;

  DefineMaterials();
  SetMaterial("Lead");
  fDetectorMessenger = new DetectorMessenger(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::~DetectorConstruction()
{ delete fDetectorMessenger;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::Construct()
{
  return ConstructVolumes();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::DefineMaterials()
{
    //Get Materials from NIST database
    G4NistManager* man = G4NistManager::Instance();
    man->SetVerbose(0);

    // define world  materials (利用NIST标准数据库)
    vacuum  = man->FindOrBuildMaterial("G4_Galactic");

  //
  // define Elements
  //
  G4double z,a;
  
  G4Element* H  = new G4Element("Hydrogen" ,"H" , z= 1., a=   1.01*g/mole);
  G4Element* N  = new G4Element("Nitrogen" ,"N" , z= 7., a=  14.01*g/mole);
  G4Element* O  = new G4Element("Oxygen"   ,"O" , z= 8., a=  16.00*g/mole);
  G4Element* Na = new G4Element("Sodium"   ,"Na", z=11., a=  22.99*g/mole);
  G4Element* Ge = new G4Element("Germanium","Ge", z=32., a=  72.59*g/mole);
  G4Element* I  = new G4Element("Iodine"   ,"I" , z=53., a= 126.90*g/mole);
  G4Element* Bi = new G4Element("Bismuth"  ,"Bi", z=83., a= 208.98*g/mole);
  
  //
  // define materials
  //
  G4double density;
  G4int ncomponents, natoms;
  G4double fractionmass;  
  
  Air =   new G4Material("Air", density= 1.290*mg/cm3, ncomponents=2);
  Air->AddElement(N, fractionmass=70.*perCent);
  Air->AddElement(O, fractionmass=30.*perCent);

  G4Material* H2l = 
  new G4Material("H2liquid", density= 70.8*mg/cm3, ncomponents=1);
  H2l->AddElement(H, fractionmass=1.);

  G4Material* H2O = 
  new G4Material("Water", density= 1.000*g/cm3, ncomponents=2);
  H2O->AddElement(H, natoms=2);
  H2O->AddElement(O, natoms=1);
  H2O->SetChemicalFormula("H_2O");
  H2O->GetIonisation()->SetMeanExcitationEnergy(78.0*eV);

  new G4Material("liquidArgon", z=18., a= 39.95*g/mole, density= 1.390*g/cm3);
  
  new G4Material("Carbon"     , z=6.,  a= 12.01*g/mole, density= 2.267*g/cm3);

  new G4Material("Aluminium"  , z=13., a= 26.98*g/mole, density= 2.700*g/cm3);

  new G4Material("Silicon"    , z=14., a= 28.09*g/mole, density= 2.330*g/cm3);
  
  new G4Material("Chromium"   , z=24., a= 51.99*g/mole, density= 7.140*g/cm3);
      
  new G4Material("Copper"     , z=29., a= 63.55*g/mole, density= 8.920*g/cm3);  

  new G4Material("Germanium"  , z=32., a= 72.61*g/mole, density= 5.323*g/cm3);

  new G4Material("Nb"         , z=41., a= 92.906*g/mole,density= 8.57*g/cm3);
    
  G4Material* NaI = 
  new G4Material("NaI", density= 3.67*g/cm3, ncomponents=2);
  NaI->AddElement(Na, natoms=1);
  NaI->AddElement(I , natoms=1);
  NaI->GetIonisation()->SetMeanExcitationEnergy(452*eV);
  
  G4Material* Iod = 
  new G4Material("Iodine", density= 4.93*g/cm3, ncomponents=1);
  Iod->AddElement(I , natoms=1);
  
  G4Material* BGO = 
  new G4Material("BGO", density= 7.10*g/cm3, ncomponents=3);
  BGO->AddElement(O , natoms=12);
  BGO->AddElement(Ge, natoms= 3);
  BGO->AddElement(Bi, natoms= 4);  

  new G4Material("Iron"       , z=26., a= 55.85*g/mole, density= 7.870*g/cm3);

  new G4Material("Silver"   ,   z=47., a=107.87*g/mole, density= 10.49*g/cm3);
  new G4Material("Stannum"   ,  z=50., a=118.71*g/mole, density= 7.37*g/cm3);
  new G4Material("Samarskite" , z=62., a=150.36*g/mole, density= 7.54*g/cm3);
  new G4Material("Tantalum" ,   z=73., a=180.95*g/mole, density= 16.69*g/cm3);
  new G4Material("Tungsten"   , z=74., a=183.85*g/mole, density= 19.25*g/cm3);
  
  new G4Material("Gold"       , z=79., a=196.97*g/mole, density= 19.30*g/cm3);  

  new G4Material("Lead"       , z=82., a=207.19*g/mole, density= 11.35*g/cm3);

  new G4Material("Uranium"    , z=92., a=238.03*g/mole, density= 18.95*g/cm3);
//  G4cout << *(G4Material::GetMaterialTable()) << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::ConstructVolumes()
{
  // Cleanup old geometry
  G4GeometryManager::GetInstance()->OpenGeometry();
  G4PhysicalVolumeStore::GetInstance()->Clean();
  G4LogicalVolumeStore::GetInstance()->Clean();
  G4SolidStore::GetInstance()->Clean();

  G4Orb * solidWorld= new G4Orb("world",worldR);
  G4LogicalVolume*      logicWorld= new G4LogicalVolume( solidWorld, vacuum, "World", 0, 0, 0);
    physiWorld = new G4PVPlacement(0,G4ThreeVector(), /*at (0,0,0)*/
                                                     logicWorld, "World",0,false,0);
  //  Must place the World Physical volume unrotated at (0,0,0).


  G4Box*   solidtarget = new G4Box("target",                        //its name
                   targetsizex/2,targetsizey/2,targetsizez/2);  //its dimensions


  logictarget = new G4LogicalVolume(solidtarget,                    //its shape
                             targetMaterial,                //its material
                             targetMaterial->GetName());    //its name

  physitarget = new G4PVPlacement(0,                         //no rotation
                           G4ThreeVector(),            //at (0,0,0)
                           logictarget,                      //its logical volume
                           targetMaterial->GetName(),       //its name
                           logicWorld,                          //its mother  volume
                           false,                      //no boolean operation
                           0);                         //copy number

  PrintParameters();


  //--------- Visualization attributes -------------------------------
  G4Colour white (1.0, 1.0, 1.0,0.2) ; // white
  logicWorld -> SetVisAttributes(new G4VisAttributes(white));
  G4Colour yellow (1.0, 1.0, 0.0) ; // yellow
  logictarget -> SetVisAttributes(new G4VisAttributes(yellow));
//    logicWorld -> SetVisAttributes(G4VisAttributes::Invisible);

  return physiWorld;
  //always return the physical World

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::PrintParameters()
{
  G4cout << "\n The Target is " << G4BestUnit(targetsizex,"Length")
         << " (L) " << G4BestUnit(targetsizex,"Length")<<" (W)"
        << G4BestUnit(targetsizez,"Length")<<" (T)"
        <<" of  "<<targetMaterial->GetName() << G4endl;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetMaterial(G4String materialChoice)
{
  // search the material by its name
  G4Material* pttoMaterial = G4Material::GetMaterial(materialChoice);
  if (pttoMaterial) {
   targetMaterial = pttoMaterial;
   UpdateGeometry();
  } else {
    G4cout << "\n--> warning from DetectorConstruction::SetMaterial : "
           << materialChoice << " not found" << G4endl;  
  } 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetSize(G4double value)
{
  targetsizez = value;
  UpdateGeometry();  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "G4RunManager.hh"

void DetectorConstruction::UpdateGeometry()
{
  if (physiWorld)
  G4RunManager::GetRunManager()->DefineWorldVolume(ConstructVolumes());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
