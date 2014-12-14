// ------------------------------------------------------------
// UCNB Neutron Beta-Decay Simulation
// Author: Brad Plaster, University of Kentucky
// ------------------------------------------------------------

#include "UCNBDetectorConstruction.hh"
#include "UCNBPhysicsList.hh"
#include "UCNBPrimaryGeneratorAction.hh"
#include "UCNBPhysicsList.hh"
#include "UCNBRunAction.hh"
#include "UCNBEventAction.hh"
#include "UCNBSteppingAction.hh"
#include "UCNBStackingAction.hh"

#include "UCNBAnalysisManager.hh"

#include "G4RunManager.hh"
#include "G4UImanager.hh"

#include "Randomize.hh"

#ifdef G4VIS_USE
#include "G4VisExecutive.hh"
#endif

#ifdef G4UI_USE
#include "G4UIExecutive.hh"
#endif

int main(int argc, char** argv)
{

  // Choose the Random engine
  //
  CLHEP::HepRandom::setTheEngine(new CLHEP::RanecuEngine);

  // Construct the default run manager
  //
  G4RunManager* runManager = new G4RunManager;

  // Creation of the analysis manager
  UCNBAnalysisManager::getInstance();

  // ROOT filename
  if (argc == 1)
    UCNBAnalysisManager::getInstance()->ROOTfilename = "UCNB.root";
  else
    UCNBAnalysisManager::getInstance()->ROOTfilename = argv[2];

  // Set mandatory initialization classes
  //
  UCNBDetectorConstruction* detector = new UCNBDetectorConstruction;
  runManager->SetUserInitialization(detector);
  //
  G4VUserPhysicsList* physics = new UCNBPhysicsList;
  runManager->SetUserInitialization(physics);

  // Set mandatory user action class
  //
  G4VUserPrimaryGeneratorAction* gen_action = new UCNBPrimaryGeneratorAction(detector);
  runManager->SetUserAction(gen_action);

  runManager->SetUserAction(new UCNBRunAction);
  runManager->SetUserAction(new UCNBEventAction);
  runManager->SetUserAction(new UCNBSteppingAction);

  G4UserStackingAction* stacking_action = new UCNBStackingAction;
  runManager->SetUserAction(stacking_action);

  // Initialize G4 kernel
  //
  runManager->Initialize();

  #ifdef G4VIS_USE
    G4VisManager* visManager = new G4VisExecutive;
    visManager->Initialize();
  #endif

  // Get the pointer to the User Interface manager
  //
  G4UImanager * UImanager = G4UImanager::GetUIpointer();

  // batch mode
  if (argc!=1)
    {
      G4String command = "/control/execute ";
      G4String fileName = argv[1];
      UImanager->ApplyCommand(command+fileName);
    }
  // interactive mode: define UI session
  else
    { 
      #ifdef G4UI_USE
        G4UIExecutive * ui = new G4UIExecutive(argc,argv);
      #ifdef G4VIS_USE
        UImanager->ApplyCommand("/control/execute vis.mac");     
      #endif
        ui->SessionStart();
        delete ui;
      #endif
    }


  // Job termination
  //
  // Free the store: user actions, physics_list and detector_description are
  //                 owned and deleted by the run manager, so they should not
  //                 be deleted in the main() program !
  //
  #ifdef G4VIS_USE
    delete visManager;
  #endif
  UCNBAnalysisManager::dispose();
  delete runManager;

  return 0;
}
