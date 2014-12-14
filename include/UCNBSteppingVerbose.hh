class UCNBSteppingVerbose;

#ifndef UCNBSteppingVerbose_h
#define UCNBSteppingVerbose_h 1

#include "G4SteppingVerbose.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class UCNBSteppingVerbose : public G4SteppingVerbose 
{
 public:
   
  UCNBSteppingVerbose();
 ~UCNBSteppingVerbose();

  void StepInfo();
  void TrackingStarted();

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
