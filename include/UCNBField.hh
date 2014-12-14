#ifndef UCNBField_h
#define UCNBField_h 1

#include "globals.hh"
#include "G4ElectroMagneticField.hh"
#include "G4FieldManager.hh"
#include "G4TransportationManager.hh"

class UCNBField
#ifndef STANDALONE
 : public G4ElectroMagneticField
#endif

{
  
public:
  UCNBField();
  void  GetFieldValue( const  double Point[4], double *Bfield ) const;
  void LoadElectricFieldMap();
  void LoadMagneticFieldMap();
	       
  G4bool DoesFieldChangeEnergy() const {return true;}

private:
  G4double zE[573];
  G4double rhoE[210];
  G4double Ez[573][210];
  G4double Erho[573][210];

  G4double zB[401];
  G4double rhoB[146];
  G4double Bz[401][146];
  G4double Brho[401][146];

};

#endif
