#include "UCNBEventAction.hh"
#include "UCNBAnalysisManager.hh"

#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4TrajectoryContainer.hh"
#include "G4Trajectory.hh"
#include "G4ios.hh"

#include "TStopwatch.h"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
 
UCNBEventAction::UCNBEventAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
 
UCNBEventAction::~UCNBEventAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
 
void UCNBEventAction::BeginOfEventAction(const G4Event*)
{
  // analysis
  UCNBAnalysisManager::getInstance()->BeginOfEvent();

  timer.Reset();
  timer.Start();
}

G4double UCNBEventAction::getEventCPUTime() {
  timer.Stop();
  G4double time = timer.CpuTime();
  timer.Continue();
  return time;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
 
void UCNBEventAction::EndOfEventAction(const G4Event* evt)
{
  timer.Stop();

  // analysis
  UCNBAnalysisManager::getInstance()->EndOfEvent();


  //G4int event_id = evt->GetEventID();
  
  // get number of stored trajectories
  //
  //G4TrajectoryContainer* trajectoryContainer = evt->GetTrajectoryContainer();
  //G4int n_trajectories = 0;
  //if (trajectoryContainer) n_trajectories = trajectoryContainer->entries();
  
  // periodic printing
  //
  //if (event_id < 100 || event_id%100 == 0) {
  //  G4cout << ">>> Event " << evt->GetEventID() << G4endl;
  //  G4cout << "    " << n_trajectories 
  //   << " trajectories stored in this event." << G4endl;
  //}
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
