#include "UCNBStackingAction.hh"

#include "UCNBAnalysisManager.hh"

#include "G4ParticleDefinition.hh"
#include "G4ParticleTypes.hh"
#include "G4Track.hh"
#include "G4ios.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

UCNBStackingAction::UCNBStackingAction()
: nCerenkovCounter(0)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

UCNBStackingAction::~UCNBStackingAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4ClassificationOfNewTrack
UCNBStackingAction::ClassifyNewTrack(const G4Track * aTrack)
{
  if(aTrack->GetDefinition() == G4OpticalPhoton::OpticalPhotonDefinition())
  { // particle is optical photon
    if(aTrack->GetParentID()>0)
    { // particle is secondary
      nCerenkovCounter++;
      G4double opticalPhotonEnergy = 0.;
      opticalPhotonEnergy = aTrack->GetKineticEnergy();
      G4double E_eV = opticalPhotonEnergy/eV;
      G4double h_eV_s = 2.*3.1416*6.582e-16;
      G4double nu_s = E_eV / h_eV_s;
      G4double wavelength = (2.998e8/nu_s) * 1.e9;
      sumCerenkovWavelength = sumCerenkovWavelength + wavelength;
    }
  }
  return fUrgent;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void UCNBStackingAction::NewStage()
{
  //G4cout << "Number of Cerenkov photons produced in this event: "
  //       << nCerenkovCounter << G4endl;

  G4double xn = (double) nCerenkovCounter;
  G4double y;
  if (xn > 0.) {
    y  = sumCerenkovWavelength / xn;
  }
  else {
    y = 0.;
  }
  //UCNBAnalysisManager::getInstance()->writeNCerPhotToNtuple(xn,y);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void UCNBStackingAction::PrepareNewEvent()
{ nCerenkovCounter = 0; sumCerenkovWavelength = 0.;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
