#ifndef UCNBSteppingAction_h
#define UCNBSteppingAction_h 1

#include "G4UserSteppingAction.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class UCNBSteppingAction : public G4UserSteppingAction
{
  public:
    UCNBSteppingAction();
    virtual ~UCNBSteppingAction(){};

    virtual void UserSteppingAction(const G4Step*);

private:

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
