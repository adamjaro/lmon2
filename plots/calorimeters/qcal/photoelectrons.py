
import ROOT as rt
from ROOT import gInterpreter

#_____________________________________________________________________________
def photoelectrons():

    #energy to wavelength PDG 2024 Tab 1.1
    gInterpreter.Declare("""
        class sens_pde {
            TGraph gp;
            public:
            void set_pde(string inp) {
                ROOT::RDataFrame df = ROOT::RDF::FromCSV(inp);
                auto df1 = df.Define("PDE_5V", "0.01*DetectionEfficiency_5V_percent");
                df1 = df1.Define("PDE_2_5V", "0.01*DetectionEfficiency_2_5V_percent");
                gp = df1.Graph("Wavelength_nm", "PDE_2_5V").GetValue();
                //gp = df1.Graph("Wavelength_nm", "PDE_5V").GetValue();
                gp.SetBit(TGraph::kIsSortedX);
            }
            vector<Bool_t> operator()(ROOT::VecOps::RVec<double>& en) {
                //cout << "sens" << endl;
                vector<Bool_t> v;
                v.resize(en.size());
                for(size_t i=0; i<en.size(); i++) {
                    v[i] = kFALSE;
                    Double_t lam = 1239.841/en[i];
                    if(lam < 300 or lam > 950) continue;
                    if( gRandom->Rndm() < 0.64*gp.Eval(lam) ) v[i] = kTRUE;
                    //cout << lam << " " << gp.Eval(lam) << " " << gRandom->Rndm() << endl;
                }
                return v;}
        };

        sens_pde sens;
    """)
    rt.sens.set_pde("../../../../lmon2-data/pde_wavelength_30035.csv")

    #number of photoelectrons in sensor in event
    gInterpreter.Declare("""
        ROOT::RVecI nphotoel_sens(const ROOT::RVecI& cell_id, const ROOT::RVecB& is_detected) {
            map<Int_t, Int_t> sens;
            for(size_t i=0; i<cell_id.size(); i++) {
                if( !is_detected[i] ) continue;
                Int_t id = cell_id[i];
                if(sens.find(id) == sens.end()) {sens.emplace(id, 1);} else {sens[id]++;}
            }
            ROOT::RVecI v;
            for(pair<Int_t, Int_t> p: sens) {v.push_back(p.second);}
            return v;
        }
    """)

    #photoelectron count in event
    gInterpreter.Declare("""
        Int_t nphotoel_evt(const ROOT::RVecB& is_detected) {
            Int_t n=0;
            for(bool i: is_detected) { if(i) n++; }
            return n;
        }
    """)

#photoelectrons


















