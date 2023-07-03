
#ifndef LoadXML_h
#define LoadXML_h

// geometry from XML

#include <string>

#include <boost/property_tree/ptree.hpp>

#include "GeoParser.h"

class LoadXML {

  public:

    LoadXML(GeoParser& geo): fGeo(geo) {}

    void ReadInput(std::string inp);

  private:

    using prop_tree = boost::property_tree::ptree;

    void ReadPTree(const prop_tree::key_type& key, const prop_tree& tree);

    void LoadDetector(const prop_tree::key_type& key, const prop_tree& tree);
    void FinishDetector();

    void LoadInclude(const prop_tree::key_type& key, const prop_tree& tree);
    void LoadIncludeDir(const prop_tree::key_type& key, const prop_tree& tree);
    void LoadConstant(const prop_tree::key_type& key, const prop_tree& tree);

    GeoParser& fGeo; // geometry to be set from loaded xml

    bool fDetectorLoad = false; // flag to indicate detector reading

    std::string fIncludeDir; // directory for include files, if set

    std::vector<std::pair<std::string, std::string>> fDetParam; // detector parameters

    std::vector<std::vector<std::string>> fDetNestedAttr; // nested attributes for detector
    std::map<std::string, int> fDetNestedAttrCount; // counter for nested attributes

    enum kdet_nes {knes_key=0, knes_cnt, knes_nam, knes_val}; // indices for nested branches

};//LoadXML

#endif

