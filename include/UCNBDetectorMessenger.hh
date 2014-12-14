#ifndef UCNBDetectorMessenger_h
#define UCNBDetectorMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"

class UCNBDetectorConstruction;
class G4UIdirectory;
class G4UIcmdWithAString;
class G4UIcmdWithADoubleAndUnit;

class UCNBDetectorMessenger: public G4UImessenger
{
  public:
    UCNBDetectorMessenger(UCNBDetectorConstruction*);
   ~UCNBDetectorMessenger();
    
    void SetNewValue(G4UIcommand*, G4String);
    
  private:
    UCNBDetectorConstruction* UCNBDetector;
    
    G4UIdirectory*             UCNBDir;
    G4UIdirectory*             measCellDir;
    G4UIcmdWithADoubleAndUnit* BFieldCmd;
};

#endif
