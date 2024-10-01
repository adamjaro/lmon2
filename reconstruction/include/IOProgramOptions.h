
#ifndef IOProgramOptions_h
#define IOProgramOptions_h

class IOProgramOptions {

  public:

    IOProgramOptions(boost::program_options::variables_map& map): fMap(map) {}

    std::unique_ptr<TChain> MakeChain(const std::string& inp, std::string tnam="DetectorTree");
    std::unique_ptr<TFile> MakeFile(const std::string& out);

    std::string GetStr(std::string par);

  private:

    boost::program_options::variables_map& fMap;

};

#endif

