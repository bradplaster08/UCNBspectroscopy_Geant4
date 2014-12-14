#ifndef UCNBDetectorConstruction_h
#define UCNBDetectorConstruction_h 1

#include "globals.hh"

#include "G4VUserDetectorConstruction.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Cons.hh"
#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"
#include "G4Material.hh"
#include "G4PVPlacement.hh"
#include "G4UserLimits.hh"

#include "G4PVParameterised.hh"

#include "G4EqMagElectricField.hh"
#include "G4PropagatorInField.hh"
#include "G4TransportationManager.hh"
#include "G4ChordFinder.hh"
#include "G4ClassicalRK4.hh"

#include "UCNBField.hh"

class UCNBDetectorConstruction : public G4VUserDetectorConstruction
{
  public:

    UCNBDetectorConstruction();
   ~UCNBDetectorConstruction();

  public:

  G4VPhysicalVolume* Construct();

  //void SetMagField(G4double);

  private:

  G4Box*             solidWorld;
  G4LogicalVolume*   logicalWorld;
  G4VPhysicalVolume* physicalWorld;

  G4Tubs*            Silicon1;
  G4LogicalVolume*   logicalSilicon1;
  G4VPhysicalVolume* physicalSilicon1;

  G4Tubs*            Silicon2;
  G4LogicalVolume*   logicalSilicon2;
  G4VPhysicalVolume* physicalSilicon2;

  G4Tubs*            Dead1;
  G4LogicalVolume*   logicalDead1;
  G4VPhysicalVolume* physicalDead1;

  G4Tubs*            Dead2;
  G4LogicalVolume*   logicalDead2;
  G4VPhysicalVolume* physicalDead2;

  G4Tubs*            Drift1;
  G4LogicalVolume*   logicalDrift1;
  G4VPhysicalVolume* physicalDrift1;

  G4Tubs*            Drift2;
  G4LogicalVolume*   logicalDrift2;
  G4VPhysicalVolume* physicalDrift2;

  G4Tubs*            DecayTrap;
  G4LogicalVolume*   logicalDecayTrap;
  G4VPhysicalVolume* physicalDecayTrap;

  G4Tubs*            DecayTrapInt;
  G4LogicalVolume*   logicalDecayTrapInt;
  G4VPhysicalVolume* physicalDecayTrapInt;


  UCNBField* Field;

  G4FieldManager *pFieldMgr;
  G4MagIntegratorStepper * pStepper;
  G4EqMagElectricField * pEquation;
  G4MagInt_Driver * pIntgrDriver;
  G4ChordFinder *pChordFinder ;
  G4PropagatorInField *propInField;

};

#endif
