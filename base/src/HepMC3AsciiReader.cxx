
//_____________________________________________________________________________
//
// generator reader for HepMC3 ascii
//
//_____________________________________________________________________________

//cmake
#include "use_hepmc3.h"
#ifdef USE_HEPMC3

//C++
#include <unordered_set>
#include <map>

//HepMC3
#include "HepMC3/ReaderAscii.h"
#include "HepMC3/GenVertex.h"
#include "HepMC3/GenParticle.h"

//Geant
#include "G4GenericMessenger.hh"
#include "G4Event.hh"
#include "G4PrimaryVertex.hh"
#include "G4PrimaryParticle.hh"
#include "G4SystemOfUnits.hh"

//local classes
#include "HepMC3AsciiReader.h"
#include "MCAttribData.h"

using namespace HepMC3;

//_____________________________________________________________________________
HepMC3AsciiReader::HepMC3AsciiReader(): G4VPrimaryGenerator(), fRead(nullptr),
  fIev(0) {

  //command for input name
  fMsg = new G4GenericMessenger(this, "/lmon/input/hepmc_ascii/");
  fMsg->DeclareProperty("name", fInputName);

  //event attributes, name in hepmc (first) and name in event output (second)
  fHepmcAttrib.insert(make_pair("truex", "true_x"));
  fHepmcAttrib.insert(make_pair("truey", "true_y"));
  fHepmcAttrib.insert(make_pair("trueQ2", "true_Q2"));
  fHepmcAttrib.insert(make_pair("trueW2", "true_W2"));
  fHepmcAttrib.insert(make_pair("true_el_pT", "true_el_pT"));
  fHepmcAttrib.insert(make_pair("true_el_theta", "true_el_theta"));
  fHepmcAttrib.insert(make_pair("true_el_phi", "true_el_phi"));
  fHepmcAttrib.insert(make_pair("true_el_E", "true_el_E"));
  fHepmcAttrib.insert(make_pair("true_el_Q2", "true_el_Q2"));
  fHepmcAttrib.insert(make_pair("Flux_[photon/s]", "flux_photon_per_s"));
  fHepmcAttrib.insert(make_pair("Power_[W]", "power_W"));
  fHepmcAttrib.insert(make_pair("num_interactions", "num_interactions"));

}//HepMC3AsciiReader

//_____________________________________________________________________________
void HepMC3AsciiReader::GeneratePrimaryVertex(G4Event *evt) {

  //open at the first call
  if(fRead == nullptr) OpenInput();

  if( (++fIev)%100000 == 0 ) {
    G4cout << "HepMC3AsciiReader::GeneratePrimaryVertex, event number: " << fIev << G4endl;
  }

  //read the event
  GenEvent mc(Units::GEV,Units::MM);
  fRead->read_event(mc);

  //event attributes
  MCAttribData *attr = new MCAttribData();

  //loop over hepmc attributes
  map<string, string>::const_iterator iatt = fHepmcAttrib.cbegin();
  for(; iatt != fHepmcAttrib.cend(); iatt++) {

    shared_ptr<DoubleAttribute> attrib = mc.attribute<DoubleAttribute>( (*iatt).first );

    //set the attribute if present
    if(attrib) {
      attr->SetVal((*iatt).second, attrib->value());
    }
  }//loop over hepmc attributes

  //set the attributes to the event
  evt->SetUserInformation(attr);

  //particles loaded from vertex loop
  unordered_set<int> used_particle_id;

  //hepmc vertex loop
  for(ConstGenVertexPtr v: mc.vertices()) {

    //create the vertex
    const FourVector vpos = v->position();
    G4PrimaryVertex *vtx = new G4PrimaryVertex(vpos.x()*mm, vpos.y()*mm, vpos.z()*mm, 0);

    //particles loop
    for(ConstGenParticlePtr p: v->particles_out()) {

      //only final outgoing particles
      if( p->status() != 1 ) continue;

      const FourVector pvec = p->momentum();

      //create the G4 particle
      G4PrimaryParticle *gp = new G4PrimaryParticle(p->pid(), pvec.px()*GeV, pvec.py()*GeV, pvec.pz()*GeV, pvec.e()*GeV);
      vtx->SetPrimary(gp);

      //mark the particle as used
      used_particle_id.insert( p->id() );

      //G4cout << "vertex particle " << p->pid() << " " << p->id() << G4endl;

    }//particles loop

    //put the vertex to the event
    evt->AddPrimaryVertex(vtx);

  }//hepmc vertex loop

  //hepmc root vertex (if present)
  G4PrimaryVertex *rvtx = 0x0;

  //loop over root vertex particles
  for(ConstGenParticlePtr p: mc.particles()) {

    //particle from root vertex
    ConstGenVertexPtr v = p->production_vertex();
    if( !v->particles_in().empty() ) continue;

    //only final outgoing particles
    if( p->status() != 1 ) continue;

    //particles not loaded from hepmc vertices
    if( used_particle_id.find( p->id() ) != used_particle_id.end() ) continue;

    //G4cout << "hepmc root particle" << " " << p->id() << G4endl;

    //create the vertex
    if( !rvtx ) {
      const FourVector vpos = v->position();
      rvtx = new G4PrimaryVertex(vpos.x()*mm, vpos.y()*mm, vpos.z()*mm, 0);
      //G4cout << vpos.x() << " " << vpos.y() << " " << vpos.z() << G4endl;
    }

    //create the G4 particle
    const FourVector pvec = p->momentum();
    G4PrimaryParticle *gp = new G4PrimaryParticle(p->pid(), pvec.px()*GeV, pvec.py()*GeV, pvec.pz()*GeV, pvec.e()*GeV);
    rvtx->SetPrimary(gp);

  }//loop over root vertex particles

  //put root vertex to the event (if present)
  if( rvtx ) {
    evt->AddPrimaryVertex(rvtx);
  }

}//GeneratePrimaryVertex

//_____________________________________________________________________________
void HepMC3AsciiReader::OpenInput() {

  G4cout << "HepMC3AsciiReader::OpenInput: " << fInputName << G4endl;

  fRead = make_shared<ReaderAscii>(fInputName.data());

}//OpenInput


#endif // cmake













