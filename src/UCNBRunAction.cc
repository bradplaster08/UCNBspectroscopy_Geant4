#include "UCNBRunAction.hh"
#include "UCNBAnalysisManager.hh"

#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4Timer.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

UCNBRunAction::UCNBRunAction()
{
  timer = new G4Timer;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

UCNBRunAction::~UCNBRunAction()
{
  delete timer;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void UCNBRunAction::BeginOfRunAction(const G4Run*)
{
  // Creation of the analysis manager
  UCNBAnalysisManager* analysis = UCNBAnalysisManager::getInstance();
  analysis->BeginOfRun();

  // Inform the runManager to save random number seed
  //G4RunManager::GetRunManager()->SetRandomNumberStore(true);
  G4RunManager::GetRunManager()->SetRandomNumberStore(false);

  // Start timer
  //G4cout << "### Run " << aRun->GetRunID() << " start." << G4endl; 
  timer->Start();

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void UCNBRunAction::EndOfRunAction(const G4Run*)
{
  // Get the analysis manager
  UCNBAnalysisManager* analysis = UCNBAnalysisManager::getInstance();
  analysis->EndOfRun();

  // Stop timer
  timer->Stop();
  G4cout << "*****G4 RUN TIMER*****" << G4endl;
  G4cout << *timer << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
