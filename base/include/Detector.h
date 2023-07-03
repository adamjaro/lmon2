
#ifndef Detector_h
#define Detector_h

// abstract base class for detectors

#include <vector>

class G4String;
class TTree;

class Detector {

  public:

    virtual ~Detector() = default;

    // add Detector to all Detectors during detector construction
    virtual void Add(std::vector<Detector*> *vec) {vec->push_back(this);}
    virtual const G4String& GetName() const = 0; // name of logical sensitive volume for Detector with sensitive volume
    virtual void CreateOutput(TTree*) {} // output of the detector

    virtual void ClearEvent() {} // beginning of event
    virtual void FinishEvent() {} // end of event

};

#endif

