#ifndef UCNBPrimaryGeneratorAction_h
#define UCNBPrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"

class UCNBDetectorConstruction;
class G4ParticleGun;
class G4Event;
 
class UCNBPrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
  public:
    UCNBPrimaryGeneratorAction(UCNBDetectorConstruction*);    
   ~UCNBPrimaryGeneratorAction();

  public:
    void GeneratePrimaries(G4Event*);

  double ENERGY2();
  double aprob2(double, double, double, double);

  private:
    G4ParticleGun* particleGun;
    UCNBDetectorConstruction* UCNBDetector;
};

#endif
