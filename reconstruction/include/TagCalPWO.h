
#ifndef TagCalPWO_h
#define TagCalPWO_h

//calorimeter analysis for CalPWO in tagger station

class TagCalPWO {

  public:

    TagCalPWO() {}

    void Run(const char *conf);

  private:

    std::string GetStr(boost::program_options::variables_map& opt_map, std::string par);



};

#endif

