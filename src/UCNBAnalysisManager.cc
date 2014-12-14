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


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "UCNBAnalysisManager.hh"
#include "G4UnitsTable.hh"

#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"

#include "UCNBStackingAction.hh"

#include "Randomize.hh"
#include "G4Poisson.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

UCNBAnalysisManager* UCNBAnalysisManager::fManager = 0;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

UCNBAnalysisManager* UCNBAnalysisManager::getInstance()
{
  if(!fManager) {
    fManager = new UCNBAnalysisManager();
  }
  return fManager;
}
void UCNBAnalysisManager::dispose()
{
  delete fManager;
  fManager = 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

UCNBAnalysisManager::UCNBAnalysisManager()
{
  verbose = 5;
  fileROOT = 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

UCNBAnalysisManager::~UCNBAnalysisManager()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void UCNBAnalysisManager::bookROOT()
{
  G4cout << "*****************" << ROOTfilename << G4endl;
  fileROOT = new TFile(ROOTfilename, "RECREATE", "UCNB Simulation");

  physics = new TTree("physics","physics");

  physics->Branch("Te0",         &Te0,         "Te0/D");
  physics->Branch("Tp0",         &Tp0,         "Tp0/D");
  physics->Branch("x0",          &x0,          "x0/D");
  physics->Branch("y0",          &y0,          "y0/D");
  physics->Branch("z0",          &z0,          "z0/D");
  physics->Branch("thetae0",     &thetae0,     "thetae0/D");
  physics->Branch("thetap0",     &thetap0,     "thetap0/D");
  physics->Branch("dEeOther",    &dEeOther,    "dEeOther/D");
  physics->Branch("dEeSilicon1", &dEeSilicon1, "dEeSilicon1/D");
  physics->Branch("dEeSilicon2", &dEeSilicon2, "dEeSilicon2/D");
  physics->Branch("dEeDead1",    &dEeDead1,    "dEeDead1/D");
  physics->Branch("dEeDead2",    &dEeDead2,    "dEeDead2/D");
  physics->Branch("dEpOther",    &dEpOther,    "dEpOther/D");
  physics->Branch("dEpSilicon1", &dEpSilicon1, "dEpSilicon1/D");
  physics->Branch("dEpSilicon2", &dEpSilicon2, "dEpSilicon2/D");
  physics->Branch("dEpDead1",    &dEpDead1,    "dEpDead1/D");
  physics->Branch("dEpDead2",    &dEpDead2,    "dEpDead2/D");
  physics->Branch("dEeSilicon1Gauss", &dEeSilicon1Gauss, "dEeSilicon1Gauss/D");
  physics->Branch("dEeSilicon2Gauss", &dEeSilicon2Gauss, "dEeSilicon2Gauss/D");
  physics->Branch("xeSilicon1",  &xeSilicon1,  "xeSilicon1/D");
  physics->Branch("yeSilicon1",  &yeSilicon1,  "yeSilicon1/D");
  physics->Branch("zeSilicon1",  &zeSilicon1,  "zeSilicon1/D");
  physics->Branch("xeSilicon2",  &xeSilicon2,  "xeSilicon2/D");
  physics->Branch("yeSilicon2",  &yeSilicon2,  "yeSilicon2/D");
  physics->Branch("zeSilicon2",  &zeSilicon2,  "zeSilicon2/D");
  physics->Branch("xpSilicon1",  &xpSilicon1,  "xpSilicon1/D");
  physics->Branch("ypSilicon1",  &ypSilicon1,  "ypSilicon1/D");
  physics->Branch("zpSilicon1",  &zpSilicon1,  "zpSilicon1/D");
  physics->Branch("xpSilicon2",  &xpSilicon2,  "xpSilicon2/D");
  physics->Branch("ypSilicon2",  &ypSilicon2,  "ypSilicon2/D");
  physics->Branch("zpSilicon2",  &zpSilicon2,  "zpSilicon2/D");
  physics->Branch("eTOF",        &eTOF,        "eTOF/D");
  physics->Branch("pTOF",        &pTOF,        "pTOF/D");
  physics->Branch("iKill",       &iKill,       "iKill/I");
  physics->Branch("timeHit1",    &timeHit1,    "timeHit1/D");
  physics->Branch("timeHit2",    &timeHit2,    "timeHit2/D");
  physics->Branch("timeDead1",   &timeDead1,   "timeDead1/D");
  physics->Branch("timeDead2",   &timeDead2,   "timeDead2/D");
  physics->Branch("eEdepTime1",  &eEdepTime1,  "eEdepTime1[500]/D");
  physics->Branch("eEdepTime2",  &eEdepTime2,  "eEdepTime2[500]/D");
  physics->Branch("eEdepDeadTime1",  &eEdepDeadTime1,  "eEdepDeadTime1[500]/D");
  physics->Branch("eEdepDeadTime2",  &eEdepDeadTime2,  "eEdepDeadTime2[500]/D");

  //physics->Branch("iStep",       &iStep,       "iStep/I");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void UCNBAnalysisManager::BeginOfRun()
{
  bookROOT();
  G4cout << ">>>>>>>>>> UCNBAnalysisManager: ROOT Ntuple is booked ..." << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void UCNBAnalysisManager::EndOfRun()
{
  saveROOT();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void UCNBAnalysisManager::BeginOfEvent()
{
  dETotal = 0.;

  dEeOther = 0.;
  dEeSilicon1 = 0.;
  dEeSilicon2 = 0.;
  dEeDead1 = 0.;
  dEeDead2 = 0.;

  dEeSilicon1Gauss = 0.;
  dEeSilicon2Gauss = 0.;

  dEpOther = 0.;
  dEpSilicon1 = 0.;
  dEpSilicon2 = 0.;
  dEpDead1 = 0.;
  dEpDead2 = 0.;

  globalTimeHit1 = -100.0*1e-9;
  globalTimeHit2 = -100.0*1e-9;

  globalTimeDead1 = -100.0*1e-9;
  globalTimeDead2 = -100.0*1e-9;

  eTOF = 0.;
  pTOF = 0.;

  xeSilicon1 = -100.;
  yeSilicon1 = -100.;
  zeSilicon1 = -100.;

  xeSilicon2 = -100.;
  yeSilicon2 = -100.;
  zeSilicon2 = -100.;

  xpSilicon1 = -100.;
  ypSilicon1 = -100.;
  zpSilicon1 = -100.;

  xpSilicon2 = -100.;
  ypSilicon2 = -100.;
  zpSilicon2 = -100.;

  iKill = -1;

  for (G4int kk=0; kk<500; kk++) {
    EdepTimeBin1[kk] = 0.;
    EdepTimeBin2[kk] = 0.;
    EdepDeadTimeBin1[kk] = 0.;
    EdepDeadTimeBin2[kk] = 0.;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void UCNBAnalysisManager::EndOfEvent()
{
  /*
  if (globalTimeHit1 > 0. && globalTimeHit2 > 0.) {
    if (globalTimeHit1 > globalTimeHit2) {
      eTOF = globalTimeHit1 - globalTimeHit2;
    }
    if (globalTimeHit2 > globalTimeHit1) {
      eTOF = globalTimeHit2 - globalTimeHit1;
    }
  }
  */
  //if (dEeSilicon1>0. && dEeSilicon2>0.) {
  //  G4cout << globalTimeHit1 << " " << globalTimeHit2 << G4endl;}
  if (globalTimeHit1 > 0. && globalTimeHit2 > 0.) {
    eTOF = globalTimeHit1 - globalTimeHit2;
  }
  else {
    eTOF = 0.;
  }


  //G4cout << "~~~~~~~~~~~~~~~~~~~~~~ " << dETotal << G4endl;

  //G4cout << globalTimeHit1 << " " << globalTimeDead1 << " " << globalTimeHit2 << " " << globalTimeDead2 << G4endl;

  timeHit1  = globalTimeHit1;
  timeHit2  = globalTimeHit2;
  timeDead1 = globalTimeDead1;
  timeDead2 = globalTimeDead2;

  for (int kk=0; kk<200; kk++) {
    //G4cout << kk << " " << EdepTimeBin1[kk] << " " << EdepTimeBin2[kk] << G4endl;
  }

  for (int kk=0; kk<500; kk++) {
    eEdepTime1[kk]     = EdepTimeBin1[kk];
    eEdepTime2[kk]     = EdepTimeBin2[kk];
    eEdepDeadTime1[kk] = EdepDeadTimeBin1[kk];
    eEdepDeadTime2[kk] = EdepDeadTimeBin2[kk];
  }


  /*
  G4double phiPotential1 = 30.;
  G4double phiPotential2 = -5.;

  if (eTOF > 0. && Te0 - ((dEeSilicon1+phiPotential1)+(dEeSilicon2)+dEeDead1+dEeDead2+dEeOther) < -1.0) {

  G4cout << "----------------------" << G4endl;
  G4cout << "PROBLEM" << G4endl;

  G4cout << Te0 << " " << cos(thetae0) << " " << eTOF/1e-9 << " " <<  dEeSilicon1 << " " << dEeDead1 << " " << dEeSilicon2 << " " << dEeDead2 << G4endl;
  G4cout << timeDead1/1e-9 << " " << timeHit1/1e-9 << " " << timeDead2/1e-9 << " " << timeHit2/1e-9 << G4endl;

  G4cout << "----------------------" << G4endl;

  }
  */

  // Electron energy resolution: assume Gaussian sigma = 2.0 keV
  // Assume: FWHM = 2.0 keV; Gaussian sigma = 2.0/2.355 = 0.85
  if (dEeSilicon1 > 0.) dEeSilicon1Gauss = G4RandGauss::shoot(dEeSilicon1, 0.85);
  if (dEeSilicon2 > 0.) dEeSilicon2Gauss = G4RandGauss::shoot(dEeSilicon2, 0.85);

  if (dEeSilicon1Gauss < 0.) dEeSilicon1Gauss = 0.;
  if (dEeSilicon2Gauss < 0.) dEeSilicon2Gauss = 0.;

  physics->Fill();
}

//void UCNBAnalysisManager::AddUpTotalEnergyDeposit(G4double x)
//{
//  dETotal = dETotal + x;
//}

void UCNBAnalysisManager::AddUpElectronOtherEnergyDeposition(G4double x)
{
  dEeOther = dEeOther + x;
}

void UCNBAnalysisManager::AddUpProtonOtherEnergyDeposition(G4double x)
{
  dEpOther = dEpOther + x;
}

void UCNBAnalysisManager::AddUpElectronSilicon1EnergyDeposition(G4double x)
{
  dEeSilicon1 = dEeSilicon1 + x;
}

void UCNBAnalysisManager::AddUpElectronSilicon2EnergyDeposition(G4double x)
{
  dEeSilicon2 = dEeSilicon2 + x;
}

void UCNBAnalysisManager::AddUpProtonSilicon1EnergyDeposition(G4double x)
{
  dEpSilicon1 = dEpSilicon1 + x;
}

void UCNBAnalysisManager::AddUpProtonSilicon2EnergyDeposition(G4double x)
{
  dEpSilicon2 = dEpSilicon2 + x;
}

void UCNBAnalysisManager::AddUpElectronDeadLayer1EnergyDeposition(G4double x)
{
  dEeDead1 = dEeDead1 + x;
}

void UCNBAnalysisManager::AddUpElectronDeadLayer2EnergyDeposition(G4double x)
{
  dEeDead2 = dEeDead2 + x;
}

void UCNBAnalysisManager::AddUpProtonDeadLayer1EnergyDeposition(G4double x)
{
  dEpDead1 = dEpDead1 + x;
}

void UCNBAnalysisManager::AddUpProtonDeadLayer2EnergyDeposition(G4double x)
{
  dEpDead2 = dEpDead2 + x;
}

void UCNBAnalysisManager::AddUpProtonDriftTime(G4double x)
{
  pTOF = pTOF + x;
}

//void UCNBAnalysisManager::AddUpElectronDriftTime(G4double x)
//{
//  eTOF = eTOF + x;
//}

void UCNBAnalysisManager::saveEventVertex(G4double Te, G4double Tp,
  G4double xvtx, G4double yvtx, G4double zvtx, G4double the, G4double thp)
{
  Te0 = Te;
  Tp0 = Tp;
  x0 = xvtx;
  y0 = yvtx;
  z0 = zvtx;
  thetae0 = the;
  thetap0 = thp;
  //G4cout << "++++++++++++++++++++++ Te0 = " << Te0 << G4endl;
}

void UCNBAnalysisManager::recordSilicon1ePosition(G4double x, G4double y,
                                                 G4double z)
{
  xeSilicon1 = x;
  yeSilicon1 = y;
  zeSilicon1 = z;  
}

void UCNBAnalysisManager::recordSilicon2ePosition(G4double x, G4double y,
                                                 G4double z)
{
  xeSilicon2 = x;
  yeSilicon2 = y;
  zeSilicon2 = z;
}

void UCNBAnalysisManager::recordSilicon1pPosition(G4double x, G4double y,
						  G4double z)
{
  xpSilicon1 = x;
  ypSilicon1 = y;
  zpSilicon1 = z;
}

void UCNBAnalysisManager::recordSilicon2pPosition(G4double x, G4double y,
						  G4double z)
{
  xpSilicon2 = x;
  ypSilicon2 = y;
  zpSilicon2 = z;
}

void UCNBAnalysisManager::killEventFlag(G4int i)
{
  iKill = i;
}

//void UCNBAnalysisManager::recordStepNumber(G4int k)
//{
//  iStep = k;
//}

void UCNBAnalysisManager::saveROOT()
{
  fileROOT->Write();
  fileROOT->Close();
}

//void UCNBAnalysisManager::writeCellEnergyDepositionToNtuple(G4double cellEnergyDeposition)
//{
//  histo->fillTuple(0, 4, cellEnergyDeposition);
//  histo->addRow(0);
//}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
