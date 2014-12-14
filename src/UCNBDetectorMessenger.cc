#include "UCNBDetectorMessenger.hh"
#include "UCNBDetectorConstruction.hh"

#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "globals.hh"

UCNBDetectorMessenger::UCNBDetectorMessenger(UCNBDetectorConstruction* myDet)
:UCNBDetector(myDet)
{ 
  UCNBDir = new G4UIdirectory("/UCNB/");
  UCNBDir->SetGuidance("UI commands specific to this example");
  
  measCellDir = new G4UIdirectory("/UCNB/MeasCell/");
  measCellDir->SetGuidance("measurement cell control");

  BFieldCmd = new G4UIcmdWithADoubleAndUnit("/UCNB/MeasCell/setBField",this);  
  BFieldCmd->SetGuidance("Set value of magnetic field");
  BFieldCmd->SetGuidance("Uniform magnetic field will be in x-direction");
  BFieldCmd->SetParameterName("Bx",false);
  BFieldCmd->SetUnitCategory("Magnetic flux density");
  BFieldCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
}

UCNBDetectorMessenger::~UCNBDetectorMessenger()
{
  delete UCNBDir;
  delete measCellDir;
  delete BFieldCmd;
}

void UCNBDetectorMessenger::SetNewValue(G4UIcommand* command,G4String newValue)
{ 
   
  //if( command == BFieldCmd )
  // { UCNBDetector->SetMagField(BFieldCmd->GetNewDoubleValue(newValue));}

}
