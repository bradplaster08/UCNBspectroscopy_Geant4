#ifndef UCNBAnalysisManager_h
#define UCNBAnalysisManager_h 1

//---------------------------------------------------------------------------
//
// ClassName:   UCNBAnalysisManager
//
// Description: Singleton class to hold analysis parameters and build histograms.
//              User cannot access to the constructor.
//              The pointer of the only existing object can be got via
//              UCNBAnalysisManager::GetInstance() static method.
//              The first invokation of this static method makes
//              the singleton object.
//
//----------------------------------------------------------------------------
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TTree.h"

#include "globals.hh"
#include <vector>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

//class UCNBHisto;

class UCNBAnalysisManager
{

public:
  // With description

  static UCNBAnalysisManager* getInstance();
  static void dispose();

private:

  UCNBAnalysisManager();
  ~UCNBAnalysisManager();

public: // Without description

  void bookROOT();
  void saveROOT();

  void BeginOfRun();
  void EndOfRun();

  void BeginOfEvent();
  void EndOfEvent();

  void saveEventVertex(G4double, G4double, G4double, G4double, G4double,
		       G4double, G4double);

  void killEventFlag(G4int);
  //void recordStepNumber(G4int);

  void AddUpElectronOtherEnergyDeposition(G4double);
  void AddUpProtonOtherEnergyDeposition(G4double);

  void AddUpElectronSilicon1EnergyDeposition(G4double);
  void AddUpElectronSilicon2EnergyDeposition(G4double);
  void AddUpElectronDeadLayer1EnergyDeposition(G4double);
  void AddUpElectronDeadLayer2EnergyDeposition(G4double);

  void AddUpProtonSilicon1EnergyDeposition(G4double);
  void AddUpProtonSilicon2EnergyDeposition(G4double);
  void AddUpProtonDeadLayer1EnergyDeposition(G4double);
  void AddUpProtonDeadLayer2EnergyDeposition(G4double);

  void recordSilicon1ePosition(G4double, G4double, G4double);
  void recordSilicon2ePosition(G4double, G4double, G4double);

  void recordSilicon1pPosition(G4double, G4double, G4double);
  void recordSilicon2pPosition(G4double, G4double, G4double);

  void AddUpTotalEnergyDeposit(G4double);
  void AddUpProtonDriftTime(G4double);

  G4int particleID;

  // electron time of flight
  G4double globalTimeHit1;
  G4double globalTimeHit2;
  G4double globalTimeDead1;
  G4double globalTimeDead2;

  // electron energy deposition vs. time
  G4int nTimeBin;
  G4double EdepTimeBin1[500], EdepTimeBin2[500];
  G4double EdepDeadTimeBin1[500], EdepDeadTimeBin2[500];

  G4String ROOTfilename;

private:

  // MEMBERS
  static UCNBAnalysisManager* fManager;

  TFile *fileROOT;
  TTree *physics;

  // ntuple tree variables
  G4double Te0, Tp0;
  G4double x0, y0, z0;
  G4double thetae0, thetap0;
  G4double dEeOther, dEeSilicon1, dEeSilicon2, dEeDead1, dEeDead2;
  G4double dEpOther, dEpSilicon1, dEpSilicon2, dEpDead1, dEpDead2;
  G4double dEeSilicon1Gauss, dEeSilicon2Gauss;
  G4double dETotal;
  G4double pTOF, eTOF;
  G4double xeSilicon1, yeSilicon1, zeSilicon1;
  G4double xeSilicon2, yeSilicon2, zeSilicon2;
  G4double xpSilicon1, ypSilicon1, zpSilicon1;
  G4double xpSilicon2, ypSilicon2, zpSilicon2;
  G4int iKill;
  //G4int iStep;
  G4double timeHit1, timeHit2, timeDead1, timeDead2;

  G4double eEdepTime1[500];
  G4double eEdepTime2[500];

  G4double eEdepDeadTime1[500];
  G4double eEdepDeadTime2[500];

  G4int verbose;


  // UCNBHisto*  histo;

  G4double dECell;
  G4double dxCell;

  G4double dEAcrylic;
  G4double dxAcrylic;

  G4double x_final, y_final, z_final;

  G4int n_step_array;
  G4double x_step_array[1000], y_step_array[1000], z_step_array[1000];
  G4double Edep_step_array[1000];

};

#endif
