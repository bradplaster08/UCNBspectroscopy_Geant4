#ifndef UCNBRunAction_h
#define UCNBRunAction_h 1

#include "G4UserRunAction.hh"
#include "globals.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class G4Run;
class G4Timer;

class UCNBRunAction : public G4UserRunAction
{
  public:
    UCNBRunAction();
   ~UCNBRunAction();

  public:
    void BeginOfRunAction(const G4Run*);
    void EndOfRunAction(const G4Run*);

  private:
    G4Timer* timer;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
