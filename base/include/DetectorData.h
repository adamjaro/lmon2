
//_____________________________________________________________________________
//
// Base class to write and read a data collection to a TTree.
//
// Objects of class Unit U holding the data are stored in Storage S container
// for writing to the tree. Derived class representing the collection
// of the Unit objects implements method to fill the Storage.
//
// The Unit objects are loaded and read from the current tree entry
// by LoadInput, GetN and GetUnit functions defined directly here.
//
//_____________________________________________________________________________

#ifndef DetectorData_h 
#define DetectorData_h

//C++
#include <vector>
#include <algorithm>
#include <string>
#include <functional>

//ROOT
#include "TTree.h"

#define DATA_ADD_UNIT_ATTR(nam) AddUnitAttr("_"#nam, fUnitIO.nam);

template<class U, class S=std::vector<U>> class DetectorData {

  // U is the Unit, S is the Storage container where objects of Unit are stored

  public:

    //_____________________________________________________________________________
    void CreateOutput(std::string base_nam, TTree *tree) {

      //make the output for the Unit storage

      //attribute representation loop
      for(UnitAttrBase *i: fUnitAttr) i->CreateOutput(base_nam, tree);

    }//CreateOutput

    //_____________________________________________________________________________
    void ClearEvent() {

      //clear run-time Unit objects and memory representation

      fStorage.clear();
      std::for_each(fUnitAttr.begin(), fUnitAttr.end(), std::mem_fn( &UnitAttrBase::ClearEvent ));

    }//ClearEvent

    //_____________________________________________________________________________
    U& Add(U obj) {

      //make new object obj of Unit U with Storage S as a sequence container
      //and return reference to it

      fStorage.push_back( obj ); // put object obj to the sequence container

      return fStorage.back(); // return reference to the added object

    }//Add

    //_____________________________________________________________________________
    template<typename key_type> U& ConstructedAt(key_type i, U obj) {

      //create object obj of unit U in the Storage S or retrieve already existing
      //with S as an associative container and return reference to it

      return ( *(fStorage.emplace(i, obj).first) ).second;

    }//ConstructedAt

    //_____________________________________________________________________________
    void FinishEvent() {

      //write all the Unit objects in the Storage to the tree

      //Unit loop
      for(auto& i: fStorage) {

        //set the IO member
        fUnitIO = GetUnitFromIter(i);

        //write the unit attributes to the tree
        std::for_each(fUnitAttr.begin(), fUnitAttr.end(), std::mem_fn( &UnitAttrBase::Write ));

      }//Unit loop

    }//FinishEvent

    //_____________________________________________________________________________
    void ConnectInput(std::string base_nam, TTree *tree) {

      //connect input for Unit attributes from the tree

      //attribute representation loop
      for(UnitAttrBase *i: fUnitAttr) i->ConnectInput(base_nam, tree);

    }//ConnectInput

    //_____________________________________________________________________________
    void LoadInput() {

      //load all Unit objects at a given tree entry

      //clear the run-time units for the event
      fUnitsRead.clear();

      //unit object loop
      for(unsigned long iu = 0; iu < fUnitAttr.front()->GetN(); iu++) {

        //attribute loop for the unit
        for(UnitAttrBase *i: fUnitAttr) {
          i->LoadVal(iu);
        }//attribute loop for the unit

        //make the unit from the loaded attributes
        fUnitsRead.push_back( U(fUnitIO) );

      }//unit object loop

    }//LoadInput

    unsigned long GetN() { return fUnitsRead.size(); } // get number of units read at the current tree entry
    U& GetUnit(unsigned long i) { return fUnitsRead[i]; } // get unit object at position i

    //read iterator
    typename std::vector<U>::iterator read_begin() { return fUnitsRead.begin(); }
    typename std::vector<U>::iterator read_end() { return fUnitsRead.end(); }

    std::vector<U>& GetReadData() { return fUnitsRead; } // handle to vector with read objects of Unit U

  protected:

    //_____________________________________________________________________________
    template<typename T> void AddUnitAttr(std::string attr_nam, T& ref) {

      //add Unit attribute of type T with its name attr_nam and reference to data value ref

      fUnitAttr.push_back( new UnitAttr<T>(attr_nam, ref) );

    }//AddUnitAttr

    //data interface for derived specific implementation
    U fUnitIO; // unit object for input/output operations
    S fStorage; // container holding the unit objects

  private:

    //unit memory representation
    class UnitAttrBase {
    public:
      virtual void CreateOutput(std::string base_nam, TTree* tree) = 0;
      virtual void ClearEvent() = 0;
      virtual void Write() = 0;
      virtual void ConnectInput(std::string base_nam, TTree* tree) = 0;
      virtual unsigned long GetN() = 0;
      virtual void LoadVal(unsigned long i) = 0;
    };

    //unit attribute in its memory representation
    template<typename T> class UnitAttr: public UnitAttrBase {
      public:
        // constructor with attribute name and value reference in fUnitIO
        UnitAttr(std::string n, T& v): attr_nam(n), val(v) {}
        void CreateOutput(std::string base_nam, TTree *tree) {
          //create the vector and tree branch from provided base name and attribute name
          vec = new std::vector<T>();
          tree->Branch((base_nam+attr_nam).c_str(), vec);
        }
        void ClearEvent() { if(vec) vec->clear(); } //clear the vector for next event
        void Write() { if(vec) vec->push_back(val); } // write the given value to the vector with attribute values
        void ConnectInput(std::string base_nam, TTree *tree) {
          //connect the vector with attribute values to the input tree
          vec = 0x0;
          if( !tree->FindBranch( (base_nam+attr_nam).c_str() ) ) return;
          tree->SetBranchAddress((base_nam+attr_nam).c_str(), &vec);
        }
        unsigned long GetN() { return vec->size(); } // number of attribute values for a given event
        void LoadVal(unsigned long i) {
          //load the value at position i via reference val
          val = 0;
          if( !vec ) return;
          val = vec->at(i);
        }
      private:
        std::string attr_nam; // attribute name
        T& val; // input/output value
        std::vector<T> *vec=0x0; // container for attribute values
    };//UnitAttr

    //Unit object from iterator for S as a sequence or an associative container
    U GetUnitFromIter(U& u_obj) {  return u_obj; } // sequence container
    template<typename some_key> U GetUnitFromIter(std::pair<some_key, U> a_obj) { return a_obj.second; } // associative container

    std::vector<UnitAttrBase*> fUnitAttr; // unit attributes in its memory representation

    std::vector<U> fUnitsRead; // unit container for reading from input tree

};//DetectorData

#endif



















