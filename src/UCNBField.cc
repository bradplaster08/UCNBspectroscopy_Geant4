#include "UCNBField.hh"
#include <iostream>
#include <fstream>
#include <TString.h>

#include <cmath>

UCNBField::UCNBField() 
{
  UCNBField::LoadElectricFieldMap();
  UCNBField::LoadMagneticFieldMap();
}

void UCNBField::LoadElectricFieldMap()
{
  TString filenameE;
  //filenameE = "/home62/phy/plaster/g4work/UCNB1.1/electric_field_z_rho.dat";
  //ifstream fileinE;
  //fileinE.open(filenameE.Data());
  //G4cout << "Reading electric field map from ... " << filenameE.Data() << G4endl;

  // Define coordinate arrays
  for (G4int i=0; i<573; i++) {
    zE[i] = ((G4double)i) - 286.;
    //G4cout << z[i] << G4endl;
  }
  for (G4int j=0; j<210; j++) {
    rhoE[j] = ((G4double)j) / 10.;
    //G4cout << rho[j] << G4endl;
  }

  // Read electric field map
  /*
  G4double ztmp, rhotmp;
  if (!fileinE) {
    G4cout << "CANNOT OPEN ELECTRIC FIELD MAP!" << G4endl;
    exit(0);
  }
  else {
    for (G4int j=0; j<210; j++) {
      for (G4int i=0; i<573; i++) {
        fileinE >> ztmp >> rhotmp >> Ez[i][j] >> Erho[i][j];
      }
    }
  }
  */
  for (G4int i=0; i<573; i++) {
    //G4cout << zE[i] << " " << Ez[i][0] << " " << Erho[i][0] << G4endl;
  }

}


void UCNBField::LoadMagneticFieldMap()
{
  TString filenameB;
  //filenameB = "/home62/phy/plaster/g4work/UCNB1.1/ami_scs_fieldmap.dat";
  //ifstream fileinB;
  //fileinB.open(filenameB.Data());
  //G4cout << "Reading magnetic field map from ... " << filenameB.Data() << G4endl;

  // Define coordinate arrays 
  for (G4int i=0; i<401; i++) {
    zB[i] = -1000. + ((G4double)i)*5.;
    //G4cout << zB[i] << G4endl;
  }
  for (G4int j=0; j<146; j++) {
    rhoB[j] = ((G4double)j)*5.;
    //G4cout << rhoB[j] << G4endl; 
  }

  // Read magnetic field map
  /*
  G4double ztmp, rhotmp, Bmag, frac;
  if (!fileinB) {
    G4cout << "CANNOT OPEN MAGNETIC FIELD MAP!" << G4endl;
    exit(0);
  }
  else {
    for (G4int j=0; j<146; j++) {
      for (G4int i=0; i<401; i++) {
        fileinB >> ztmp >> rhotmp >> Bz[i][j] >> Brho[i][j] >> Bmag >> frac;
      }
    }
  }
  */
  for (G4int i=0; i<401; i++) {
    //G4cout << zB[i] << " " << Bz[i][0] << " " << Brho[i][0] << G4endl;
  }

}

void UCNBField::GetFieldValue(const double point[4], double *Bfield ) const
{
  // These (x,y,z) positions are returned in units of [mm]
  G4double x_mm = point[0];
  G4double y_mm = point[1];
  G4double z_mm = point[2];

  // Convert to units of [m]
  G4double x_cm, y_cm, z_cm;
  G4double x_m = x_mm/1000.;
  G4double y_m = y_mm/1000.;
  G4double z_m = z_mm/1000.;

  // Sinusoidal interpolation of field
  G4double zgrid[6];
  zgrid[0] = -9.0;
  zgrid[1] = -2.2;
  zgrid[2] = -1.5;
  zgrid[3] =  1.5;
  zgrid[4] =  2.2;
  zgrid[5] =  9.0;

  G4double Bzgrid[6];
  Bzgrid[0] = 0.6;
  Bzgrid[1] = 0.6;
  Bzgrid[2] = 1.0;
  Bzgrid[3] = 1.0;
  Bzgrid[4] = 0.6;
  Bzgrid[5] = 0.6;

  int zindexLow = 0;
  for (G4int i=0; i<5; i++) {
    if ( (z_m > zgrid[i]) && (z_m < zgrid[i+1]) ) zindexLow = i;
  }

  //G4cout << z_m << " " << zindexLow << G4endl;

  G4double base = (Bzgrid[zindexLow] + Bzgrid[zindexLow+1]) / 2.;
  G4double amp  = (Bzgrid[zindexLow] - Bzgrid[zindexLow+1]) / 2.;
  G4double dz   = zgrid[zindexLow+1] - zgrid[zindexLow];
  G4double l    = (z_m - zgrid[zindexLow])/dz;

  /*
  G4double B0Tesla;
  if ((z_m>2.) || (z_m<-2.)) B0Tesla = 1.0;
  else B0Tesla = 1.0; // B0 in Tesla

  G4double E0; // E0 field in kV/cm
  if (z_m>=2.00 && z_m<=2.10) E0 = 3.0;
  else if (z_m>=-2.10 && z_m<=-2.00) E0 = -3.0;
  else E0 = 0.0;
  */

  G4double Bxint, Byint, Bzint, Brhoint;
  G4double Exint, Eyint, Ezint, Erhoint;

  Bzint = base + amp*cos(l*M_PI);

  Brhoint = amp*M_PI*sin(l*M_PI)/(2.*dz);
  Bxint   = Brhoint * x_m;
  Byint   = Brhoint * y_m;

  // Magnetic field: (x,y,z) -> (0,1,2) components of Bfield[.]
  Bfield[0] = Bxint * tesla;
  Bfield[1] = Byint * tesla;
  Bfield[2] = Bzint * tesla;
  //Bfield[0] = 0. * tesla;
  //Bfield[1] = 0. * tesla;
  //Bfield[2] = 0. * tesla;

  /*
  if (z_m > -2.2 && z_m < -1.5) {
  G4cout << x_m << " " << y_m << " " << z_m << " " << Bfield[0]/tesla << " " << Bfield[1]/tesla << " " << Bfield[2]/tesla << G4endl;
  }
  */

  G4double deltaphiNeg =  -5.; // [kV]
  G4double deltazNeg   =  20.; // [cm]

  G4double deltaphiPos = -30.; // [kV]
  G4double deltazPos   =  20.; // [cm]

  Exint = 0.;
  Eyint = 0.;
  Ezint = 0.;

  G4double tiny = 1.0e-9;

  // asymmetric_run1:
  //if ( (z_m > -2.2 + tiny) && (z_m <= -2.0) ) Ezint = -1.*deltaphiNeg/deltazNeg;
  //if ( (z_m >=  2.0) && (z_m < 2.2 - tiny) ) Ezint = -1.*deltaphiPos/deltazPos;

  // asymmetric_run2:
  if ( (z_m >= -1.5) && (z_m <= 1.5) ) Ezint = 0.1;
  //if ( (z_m >= -1.5) && (z_m <= 1.5) ) Ezint = 0.15;
  if ( (z_m > 1.5) && (z_m < 2.2 - tiny) ) Ezint = (30./70.);

  // Electric field: (x,y,z) -> (3,4,5) components of Bfield[.]
  Bfield[3] = Exint * kilovolt/cm;
  Bfield[4] = Eyint * kilovolt/cm;
  Bfield[5] = Ezint * kilovolt/cm;

  //G4cout << Bfield[0] << " " << Bfield[3] << G4endl;

}
