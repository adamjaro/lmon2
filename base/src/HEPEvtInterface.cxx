
//_____________________________________________________________________________
//
// reader for HEPEvt, based on G4HEPEvtInterfaceIons
//_____________________________________________________________________________

//C++ headers
#include <string.h>

//Geant headers
#include "G4GenericMessenger.hh"
#include "G4SystemOfUnits.hh"
#include "G4PrimaryVertex.hh"
#include "G4Event.hh"
#include "G4IonTable.hh"

//local headers
#include "HEPEvtInterface.h"

using namespace std;

//_____________________________________________________________________________
HEPEvtInterface::HEPEvtInterface() : G4VPrimaryGenerator(), fInputName(""), vLevel(0), fIev(0) {

  //command for name of input file
  fMsg = new G4GenericMessenger(this, "/lmon/input/hepevt/");
  fMsg->DeclareProperty("name", fInputName);
  fMsg->DeclareProperty("vLevel", vLevel);

}//HEPEvtInterface

//_____________________________________________________________________________
void HEPEvtInterface::GeneratePrimaryVertex(G4Event *evt) {

  //open input at the first call
  if(!fIn.is_open()) OpenInput();

  //increment event count for progress printout
  fIev++;
  if( fIev%100000 == 0 ) {
    G4cout << "HEPEvtInterface::GeneratePrimaries, event number: " << fIev << G4endl;
  }

  //number of entries
  G4int NHEP = 0;
  fIn >> NHEP;

  if( !fIn.good() ) {
    G4cout << "HEPEvtInterface::GeneratePrimaries: no more events" << G4endl;
    return;
  }

  //particle loop
  for( G4int IHEP=0; IHEP<NHEP; ++IHEP )
  {
    G4int ISTHEP;   // status code
    G4int IDHEP;    // PDG code
    G4int JDAHEP1;  // first daughter
    G4int JDAHEP2;  // last daughter
    G4double PHEP1; // px in GeV
    G4double PHEP2; // py in GeV
    G4double PHEP3; // pz in GeV
    G4double PHEP5; // mass in GeV

    fIn >> ISTHEP >> IDHEP >> JDAHEP1 >> JDAHEP2
              >> PHEP1 >> PHEP2 >> PHEP3 >> PHEP5;
    if( fIn.eof() ) {
      G4Exception("HEPEvtInterface::GeneratePrimaryVertex", "Event0203",
                  FatalException,
                  "Unexpected End-Of-File in the middle of an event");
    }
    if(vLevel > 1)
    {
      G4cout << " " << ISTHEP << " " << IDHEP << " " << JDAHEP1
             << " " << JDAHEP2 << " " << PHEP1 << " " << PHEP2
             << " " << PHEP3 << " " << PHEP5 << G4endl;
    }

    ///Create G4PrimaryParticle object
    G4PrimaryParticle* particle = nullptr;
    if( IDHEP < 1000000000 ) {
      //PDG particle

      particle = new G4PrimaryParticle( IDHEP );
    } else {
      //ion

      std::string ipdg = std::to_string(IDHEP);
      int atomic_number = std::stoi( ipdg.substr(3,3) );
      int atomic_mass = std::stoi( ipdg.substr(6,3) );
      int excit_level = std::stoi( ipdg.substr(9,1) );
      G4ParticleDefinition *ion = G4IonTable::GetIonTable()->GetIon(atomic_number, atomic_mass, excit_level);
      particle = new G4PrimaryParticle(ion);
    }
    particle->SetMass( PHEP5*GeV );
    particle->SetMomentum(PHEP1*GeV, PHEP2*GeV, PHEP3*GeV );

    // Create G4HEPEvtParticle object
    //
    G4HEPEvtParticle* hepParticle = new G4HEPEvtParticle( particle, ISTHEP, JDAHEP1, JDAHEP2 );

    // Store
    fHPlist.push_back( hepParticle );
  }

  //at least one particle
  if( fHPlist.size() == 0 ) return; 

  //connection between daughter particles decayed from the same mother
  for( std::size_t i=0; i<fHPlist.size(); ++i )
  {
    if( fHPlist[i]->GetJDAHEP1() > 0 ) //  it has daughters
    {
      G4int jda1 = fHPlist[i]->GetJDAHEP1()-1; // FORTRAN index starts from 1
      G4int jda2 = fHPlist[i]->GetJDAHEP2()-1; // but C++ starts from 0.
      G4PrimaryParticle* mother = fHPlist[i]->GetTheParticle();
      for( G4int j=jda1; j<=jda2; ++j )
      {
        G4PrimaryParticle* daughter = fHPlist[j]->GetTheParticle();
        if(fHPlist[j]->GetISTHEP()>0)
        {
          mother->SetDaughter( daughter );
          fHPlist[j]->Done();
        }
      }
    }
  }

  // Create G4PrimaryVertex object
  //
  G4PrimaryVertex* vertex = new G4PrimaryVertex(0., 0., 0., 0.);

  //initial particles in vertex
  for( std::size_t ii=0; ii<fHPlist.size(); ++ii )
  {
    if( fHPlist[ii]->GetISTHEP() > 0 ) // ISTHEP of daughters had been 
                                      // set to negative
    {
      G4PrimaryParticle* initialParticle = fHPlist[ii]->GetTheParticle();
      vertex->SetPrimary( initialParticle );
    }
  }

  // Clear G4HEPEvtParticles
  //
  for(std::size_t iii=0; iii<fHPlist.size(); ++iii)
  {
    delete fHPlist[iii];
  }
  fHPlist.clear();

  // Put the vertex to G4Event object
  //
  evt->AddPrimaryVertex( vertex );

}//GeneratePrimaryVertex

//_____________________________________________________________________________
void HEPEvtInterface::OpenInput() {

  //open the input file

  G4cout << "HEPEvtInterface::OpenInput: " << fInputName << G4endl;
  fIn.open(fInputName);

  //test if file exists
  if(fIn.fail()) {
    string description = "Can't open input: " + fInputName;
    G4Exception("HEPEvtInterface::OpenInput", "InputNotOpen01", FatalException, description.c_str());
  }

}//OpenInput























