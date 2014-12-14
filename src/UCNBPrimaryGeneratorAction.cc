#include "UCNBPrimaryGeneratorAction.hh"
#include "UCNBDetectorConstruction.hh"
#include "UCNBAnalysisManager.hh"

#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "globals.hh"
#include "Randomize.hh"

UCNBPrimaryGeneratorAction::UCNBPrimaryGeneratorAction(
                                                       UCNBDetectorConstruction* myDC)
:UCNBDetector(myDC)
{
}

UCNBPrimaryGeneratorAction::~UCNBPrimaryGeneratorAction()
{
  delete particleGun;
}

void UCNBPrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{

  // Robby's electron-proton event generator

  G4double MN = 939.565379*1e6; // eV
  G4double MP = 938.272046*1e6; // eV
  G4double ME = 510.998928*1e3; // eV
  G4double elambda = -1.2701;
  G4double alpha = 1./137.036;

  G4double POL = 1.0;
  G4double reject = 100.0;

  G4double MN1 = MN/ME;
  G4double MP1 = MP/ME;
  G4double ME1 = 1.;
  G4double DELTA = MN1 - MP1;

  G4double XCRIT = 0.5*(DELTA + ME1*ME1/DELTA);
  G4double EEMAX = DELTA - (DELTA*DELTA-1.)/(2.*MN1);

  G4double A = -2.*elambda*(1.+elambda)/(1.+3.*elambda*elambda);
  G4double B = -2.*elambda*(1.-elambda)/(1.+3.*elambda*elambda);
  G4double ALIT = (1.-elambda*elambda)/(1.+3.*elambda*elambda);
  G4double el1 = elambda;
  G4double BLIT = 0.0;

  G4double gg = 1.08;
  G4double ff = 0.00;
  G4double EE, PE, NEPHI, NTHETA, NEU, NEV, NEW, THETAE, PHIE;

  //G4cout << " " << G4endl;
  //G4cout << "Starting..." << G4endl;

  while (gg>ff) {
    EE = ENERGY2() + 1.;
    PE = sqrt(EE*EE-1.);

    NEPHI = 0.0;
    NTHETA = 0.0;

    NEU = cos(NEPHI)*sin(NTHETA);
    NEV = sin(NEPHI)*sin(NTHETA);
    NEW = cos(NTHETA);

    THETAE = M_PI*G4UniformRand();
    PHIE = 2.*M_PI*G4UniformRand();

    ff = aprob2(A,el1,THETAE,EE);
    gg = 1.08*G4UniformRand();
  }

  //G4cout << "Te = " << (EE-1.)*ME/1000. << G4endl;

  G4double EU = cos(PHIE)*sin(THETAE);
  G4double EV = sin(PHIE)*sin(THETAE);
  G4double EW = cos(THETAE);

  G4double THETAEP = 1.;
  G4double THETAMAX = 0.;
  G4double THETAP, PHIP, PU, PV, PW;

  tryAgain:
  THETAP = acos(2.*G4UniformRand()-1);

  if ( std::abs(THETAP - M_PI/2.) < 0.0001 ) {
    if (THETAP > M_PI/2.) {THETAP = THETAP + 0.001;}
    else {THETAP = THETAP - 0.001;}
  }

  PHIP = 2.*M_PI*G4UniformRand();
  PU = cos(PHIP)*sin(THETAP);
  PV = sin(PHIP)*sin(THETAP);
  PW = cos(THETAP);
  THETAEP = EU*PU+PV*EV+EW*PW;
  if (EE >= XCRIT) {
    THETAMAX = -1.*sqrt(1.-((MP1/MN1)*(EEMAX-EE)/PE)*((MP1/MN1)*(EEMAX-EE)/PE));
    if (THETAEP>THETAMAX) goto tryAgain;
  }

  G4double XM = MN1 - EE;
  G4double X = XM*XM-PE*PE+MP1*MP1;
  G4double Y = PE*THETAEP;

  G4double aa = 4.*(XM*XM-Y*Y);
  G4double bb = 4.*Y*X;
  G4double hbig = 0.5*(MP1+MN1+1./(MP1+MN1));
  G4double hlit = 4.*(MN1*MN1-MP1*MP1)*(EE-XCRIT)*(hbig-EE);
  G4double s = sqrt(bb*bb - 4.*aa*hlit);

  G4double PPP = (-bb + s)/(2.*aa);
  G4double PPM = (-bb - s)/(2.*aa);

  G4double RRP = sqrt(PPP*PPP+MP1*MP1);
  G4double RRM = sqrt(PPM*PPM+MP1*MP1);

  G4double EPMAX = sqrt((PE+(MN1*(DELTA-EE)/(MN1-EE-PE)))*
			(PE+(MN1*(DELTA-EE)/(MN1-EE-PE))) + MP1*MP1);

  G4double EPMIN;
  if (EE < XCRIT) {
    EPMIN = sqrt(((MN1*(DELTA-EE)/(MN1-EE-PE))-PE)*
                 ((MN1*(DELTA-EE)/(MN1-EE-PE))-PE) + MP1*MP1);
  }
  else {
    EPMIN = sqrt((PE-(MN1*(DELTA-EE)/(MN1-EE-PE)))*
                 (PE-(MN1*(DELTA-EE)/(MN1-EE-PE))) + MP1*MP1);
  }

  G4int jblah2, jblah3;
  if (EE>XCRIT) {
    if (RRP<EPMAX && RRP>EPMIN) jblah2 = 1;
    if (RRM<EPMAX && RRM>EPMIN) jblah3 = 1;
  }

  G4double FDP = std::abs(PPP/RRP + (PPP+Y)/(MN1-EE-RRP));
  G4double FDM = std::abs(PPM/RRM + (PPM+Y)/(MN1-EE-RRM));

  G4double E1P = sqrt(PPP*PPP + PE*PE + 2.*PE*PPP*THETAEP);
  G4double E1M = sqrt(PPM*PPM + PE*PE + 2.*PE*PPM*THETAEP);

  G4double CHI1 =  ( 1. - ALIT*(PE*PE+Y*PPP)/(EE*E1P) + BLIT/EE 
		     + POL*(A*PE/EE*EW-B*(PE*EW+PPP*PW)/E1P))
    * EE*PE*PPP*PPP*sin(THETAP)*sin(THETAE);

  G4double CHI2 = ( 1. - ALIT*(PE*PE+Y*PPM)/(EE*E1M) + BLIT/EE
		    + POL*(A*PE/EE*EW-B*(PE*EW+PPM*PW)/E1M))
    * EE*PE*PPM*PPM*sin(THETAP)*sin(THETAE);

  G4double DIFFP = 1.;
  G4double DIFFM = 1.;

  G4double H, F1, F2, EP, PP;
  if (EE >= XCRIT) {
    H = G4UniformRand();
    F1 = (CHI1*DIFFP/FDP)/(CHI1*DIFFP/FDP + CHI2*DIFFM/FDM);
    F2 = (CHI2*DIFFM/FDM)/(CHI1*DIFFP/FDP + CHI2*DIFFM/FDM);
    if (H > F2) {EP = RRP; PP = PPP;}
    else {EP = RRM; PP = PPM;}
  }
  else if (EE < XCRIT) {
    EP = RRP; PP = PPP; CHI2 = 0.;
  }
  if ( (EP > EPMAX) || (EP<EPMIN) ) goto tryAgain;

  G4double BETA = PE/EE;
  G4double BETAR = std::abs(BETA-(1.-BETA*BETA)*PP/EP*THETAEP);
  G4double gam = 0.5772;

  G4double FERMI = 1. + alpha*M_PI/BETAR + alpha*alpha*
    (11./4.- gam - log(2.*BETAR*EE*(0.01)/4.) + M_PI*M_PI/(3.*BETAR*BETAR));

  G4double DGAMMA1 = FERMI*(CHI1*DIFFP/FDP+CHI2*DIFFM/FDM);

  G4double REJECT;
  if (EE > XCRIT) {REJECT = 700.;}
  else if (EE <= XCRIT) {REJECT = 10.;}
  G4double G = REJECT*G4UniformRand();

  if (G > DGAMMA1) goto tryAgain;

  G4double EN = MN1 - EE - EP;
  G4double NU = -1.*(PE*EU+PP*PU)/EN;
  G4double NV = -1.*(PE*EV+PP*PV)/EN;
  G4double NW = -1.*(PE*EW+PP*PW)/EN;
  G4double THETANE = EU*NU+EV*NV+EW*NW;

  // Electron variables
  G4double Te0 = (((EE - 1.) * ME)/1000.)*keV;
  G4double px_hat_e = EU;
  G4double py_hat_e = EV;
  G4double pz_hat_e = EW;

  // Proton variables
  G4double Tp0 = (EP*ME - MP)*eV;
  G4double px_hat_p = PU;
  G4double py_hat_p = PV;
  G4double pz_hat_p = PW;

  //G4double testtheta = M_PI/2 + 0.0001;
  //G4double testphi = 2*M_PI*G4UniformRand();
  //px_hat_p = sin(testtheta)*cos(testphi);
  //py_hat_p = sin(testtheta)*sin(testphi);
  //pz_hat_p = cos(testtheta);

  // Sample vertex position
  G4double x_test, y_test;
  G4double x_vertex, y_vertex, z_vertex;

  G4double testRadius = 1.e99;
  while (testRadius > 0.062) {
    x_test = (-0.062 + G4UniformRand()*2.*0.062)*m;
    y_test = (-0.062 + G4UniformRand()*2.*0.062)*m;
    testRadius = sqrt(x_test*x_test + y_test*y_test);
  }
  x_vertex = x_test*m;
  y_vertex = y_test*m;
  //z_vertex = (-1.5 + G4UniformRand()*3.0)*m;
  z_vertex = (-1.4 + G4UniformRand()*2.8)*m;

  //x_vertex = 0.0*m;
  //y_vertex = 0.0*m;
  //z_vertex = 0.0*m;

  //G4double testAngle = acos(pz_hat_e);
  //Te0 = 434.5*keV;
  //pz_hat_e = 0.28;
  //Tp0 = 500.*eV;

  // Generate electron and proton
  G4int n_particle = 1;
  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();

  particleGun = new G4ParticleGun(n_particle);
  G4ParticleDefinition* particle1 = particleTable->FindParticle("e-");
  particleGun->SetParticleDefinition(particle1);
  //px_hat_e = 0.; py_hat_e = 0.; pz_hat_e = 1.0;
  particleGun->SetParticleMomentumDirection(G4ThreeVector(px_hat_e,py_hat_e,pz_hat_e));
  particleGun->SetParticleEnergy(Te0);
  particleGun->SetParticlePosition(G4ThreeVector(x_vertex,y_vertex,z_vertex));
  particleGun->GeneratePrimaryVertex(anEvent);

  particleGun = new G4ParticleGun(n_particle);
  G4ParticleDefinition* particle2 = particleTable->FindParticle("proton");
  particleGun->SetParticleDefinition(particle2);
  //px_hat_p = 0.; py_hat_p = 0.; pz_hat_p = -1.0;
  particleGun->SetParticleMomentumDirection(G4ThreeVector(px_hat_p,py_hat_p,pz_hat_p));
  particleGun->SetParticleEnergy(Tp0);
  particleGun->SetParticlePosition(G4ThreeVector(x_vertex,y_vertex,z_vertex));
  particleGun->GeneratePrimaryVertex(anEvent);

  G4double thetaElectron = acos(pz_hat_e);
  G4double thetaProton   = acos(pz_hat_p);

  // Save initial vertex variables
  UCNBAnalysisManager::getInstance()->saveEventVertex(Te0/keV, Tp0/keV,
    x_vertex/m, y_vertex/m, z_vertex/m, thetaElectron, thetaProton);

  //G4cout << "***ANGLES***" << G4endl;
  //G4cout << thetaProton << G4endl;
  //G4cout << thetaElectron << G4endl;
}




double UCNBPrimaryGeneratorAction::ENERGY2(){
  while( 1 ) {
    //G4double E0 = 2.5309874e0;
    G4double E0 = 2.5295196;
    G4double b = 0.0;

    G4double f = 0.00;
    G4double y = 1.80;

    G4double E, GAM, FERMI;
    while (f<y) {
      E=(E0-1.0)*G4UniformRand();
      y=1.80*G4UniformRand();
      GAM=(-E/(sqrt(E*E+2*E)*137.));  // numerator should be E+1
      //FERMI=2.*3.14*GAM/(exp(2.*3.14*GAM)-1.);
      FERMI = 1.;
      f=FERMI*sqrt(E*E+2*E)*(E0-(E+1))*(E0-(E+1))*(E+1)*(1+b*1/(E+1));
    }

    // T_e = 150, 300, 600, 650 keV
    //G4double fixedT = 750.;
    //E = fixedT / 510.998928;

    return E;

  }
}

double UCNBPrimaryGeneratorAction::aprob2(double A, double lambda,
					  double theta, double E) {
  while( 1 ) {
    double beta = sqrt(E*E - 1.)/E;
    double LAMBDA = -1.0*lambda;
    double W = 2.529476;
    double alpha = 1./137.036;
    double MM = 1837.0;
    double MU = 3.71;

    double AM = (LAMBDA+MU)/(LAMBDA*MM*(1.0-LAMBDA)*(1.0+3.*LAMBDA*LAMBDA));
    double A1 = LAMBDA*LAMBDA+(2./3.)*LAMBDA-(1./3.);
    double A2 = -1.*LAMBDA*LAMBDA*LAMBDA-3.*LAMBDA*LAMBDA-
      (5./3.)*LAMBDA+(1./3.);
    double A3 = (1.0-LAMBDA)*2.*LAMBDA*LAMBDA;

    //double Abb = A*(1.0+AM*(A1*W+A2*E+A3/E));
    double Abb = A;
    double temp = sin(theta)*(1.0+Abb*beta*cos(theta));

    return temp;
  }
}
