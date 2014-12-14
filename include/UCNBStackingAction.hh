#ifndef UCNBStackingAction_H
#define UCNBStackingAction_H 1

#include "globals.hh"
#include "G4UserStackingAction.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class UCNBStackingAction : public G4UserStackingAction
{
  public:
    UCNBStackingAction();
   ~UCNBStackingAction();

  public:
    G4ClassificationOfNewTrack ClassifyNewTrack(const G4Track* aTrack);
    void NewStage();
    void PrepareNewEvent();

  private:
    G4int nCerenkovCounter;
    G4double sumCerenkovWavelength;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

