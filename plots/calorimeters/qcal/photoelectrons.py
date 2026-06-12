
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
                ROOT::RDataFrame df = ROOT::RDF::FromCSV(inp); // input PDE vs. wavelength
                auto df1 = df.Define("PDE_5V", "0.01*DetectionEfficiency_5V_percent"); // convert from percent
                df1 = df1.Define("PDE_2_5V", "0.01*DetectionEfficiency_2_5V_percent");
                gp = df1.Graph("Wavelength_nm", "PDE_2_5V").GetValue(); // PDE as TGraph
                //gp = df1.Graph("Wavelength_nm", "PDE_5V").GetValue();
                gp.SetBit(TGraph::kIsSortedX); // sorted flag for faster evaluation
            }
            //detected flags for photons in event by energy
            vector<Bool_t> operator()(ROOT::VecOps::RVec<double>& en) {
                //cout << "sens" << endl;
                vector<Bool_t> v; // flags to return
                v.resize(en.size());
                for(size_t i=0; i<en.size(); i++) {
                    v[i] = kFALSE; // initialize as not detected
                    Double_t lam = 1239.841/en[i]; // energy (eV) to wavelength (nm), PDG 2024 Tab 1.1
                    if(lam < 300 or lam > 950) continue; // range (nm) for input PDE

                    //uniform random [0,1] compared to PDE at a given wavelength, detected flag when within
                    if( gRandom->Rndm() < gp.Eval(lam) ) v[i] = kTRUE;
                    //cout << lam << " " << gp.Eval(lam) << " " << gRandom->Rndm() << endl;
                }
                //export the detected flags
                return v;}
        };

        sens_pde sens; // PDE callable object
    """)
    rt.sens.set_pde("../../../../lmon2-data/pde_wavelength_30035.csv") # set the input PDE

    #number of photoelectrons in sensor in event
    gInterpreter.Declare("""
        ROOT::RVecI nphotoel_sens(const ROOT::RVecI& cell_id, const ROOT::RVecB& is_detected) {
            map<Int_t, Int_t> sens; // sensor ID and photoelectron count
            for(size_t i=0; i<cell_id.size(); i++) {
                if( !is_detected[i] ) continue; // test for detected
                Int_t id = cell_id[i]; // sensor ID

                //increment count for the given sensor
                if(sens.find(id) == sens.end()) {sens.emplace(id, 1);} else {sens[id]++;}
            }

            //export the counts per sensor
            ROOT::RVecI v;
            for(pair<Int_t, Int_t> p: sens) {v.push_back(p.second);} // loop over sensor, take individual counts
            return v;
        }
    """)

    #photoelectron count in event
    gInterpreter.Declare("""
        Int_t nphotoel_evt(const ROOT::RVecB& is_detected) {
            Int_t n=0;
            for(bool i: is_detected) { if(i) n++; } // sum over all detected
            return n;
        }
    """)

    #sensor count with photoelectrons above or equal a threshold
    gInterpreter.Declare("""
        Int_t nsens_photoel(const ROOT::RVecI& nphotoel, Int_t thres) {
            Int_t n=0; // sensor count to return
            for(Int_t i: nphotoel) { if(i >= thres) n++; } // increment sensor count when at least at threshold
            return n;
        }
    """)

    #photoelectron time on individual sensors in event
    gInterpreter.Declare("""
        class sensor_time {
            map<Int_t, vector<double>> sens; // sensor ID and photoelectrons time
            vector<Int_t> sens_id;
            public:
            int operator()(const ROOT::RVecI& cell_id, const ROOT::RVecD& tim, const ROOT::RVecB& is_detected) {
                //cout << "hi from sensor_time" << endl;
                for(size_t i=0; i<is_detected.size(); i++) {
                    if( !is_detected[i] ) continue; // detected flag
                    //cout << cell_id[i] << " " << tim[i] << endl;
                    sens[cell_id[i]].push_back(tim[i]);
                }
                for(const pair<Int_t, vector<double>>&i: sens) {
                    //cout << i.first << " " << i.second.size() << endl;
                    sens_id.push_back(i.first);
                }
                //cout << sens.size() << endl;
                return 0;
            }
            vector<Int_t>& GetSensId() { return sens_id; }
            map<Int_t, vector<double>>& GetSens() { return sens; }
        };
        sensor_time senstime;
    """)

    #cell position in x and z
    gInterpreter.Declare("""
        class cell_pos_xz {
            public:
            map<Int_t, pair<Double_t, Double_t>> xzpos;
            cell_pos_xz(Int_t nx, Int_t nz, Double_t cell_xy, Double_t cell_z, Double_t cell_phi, Double_t modx, Double_t modz) {
                cout << "hi from cell_pos_xz " << nx << " " << nz << " " << cell_phi << " " << modz << endl;
                Double_t cell_posz_init = 0.5*(cell_xy*TMath::Cos(cell_phi) + cell_z*TMath::Sin(cell_phi));
                Double_t cell_posz = 0.5*modz - cell_posz_init;
                Int_t cell_cnt = 0;
                //z-loop
                for(Int_t iz=0; iz<nz; iz++) {
                    //x-loop
                    for(Int_t ix=0; ix<nx; ix++) {

                        Double_t cell_posx = -0.5*modx + 0.5*cell_xy + ix*cell_xy;

                        xzpos.emplace(cell_cnt++, make_pair(cell_posx, cell_posz));
                    }//x-loop
                    cell_posz -= cell_xy/TMath::Cos(cell_phi);
                }//z-loop
                //for(size_t i=0; i<xzpos.size(); i++) cout << i << " " << xzpos[i].first << " " << xzpos[i].second << endl;
            }
            ROOT::RVecD operator()(const ROOT::RVecI& cell_id, const ROOT::RVecB& is_detected, Int_t xz=0) {
                //vector to return, 0: x, 1: z
                ROOT::RVecD vec;
                for(size_t i=0; i<cell_id.size(); i++) {
                    if( !is_detected[i] ) continue;
                    //cout << cell_id[i] << endl;
                    if( xz == 0) {
                        vec.push_back( xzpos[cell_id[i]].first );
                    } else {
                        vec.push_back( xzpos[cell_id[i]].second );
                    }
                }
                return vec;
            }
        };
    """)


#photoelectrons


















