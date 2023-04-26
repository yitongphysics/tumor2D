// header files

#include "tumor3D.h"
#include <sstream>

// preprocessor macros
#define NDIM 3

// namspace
using namespace std;

const double gamtt = 0.0;                 // surface tension
const double boxLengthScale = 2.0;        // neighbor list box size in units of initial l0
const double dt0 = 0.1;                // initial magnitude of time step in units of MD time
const double Ftol = 1e-7;

int main(int argc, char const *argv[])
{
    // forces for simulations
    tumor3DMemFn invasionForceUpdate = nullptr;
    invasionForceUpdate = &tumor3D::repulsiveTumorInterfaceForceUpdate;

    // local variables to be read in
    int NT, NPRINTSKIP, seed;
    double NTdbl, NPRINTSKIPdbl, l1, l2, v0, Dr0, aCalA0, kecm, ecmbreak, dDr, dPsi, Drmin, kv, ka, kb, kc, M, P0, volumeRatio;

    
    // read in parameters from command line input
    string inputFile         = argv[1];                // input file with initial configuration
    string NT_str             = argv[2];                // # of time steps
    string NPRINTSKIP_str     = argv[3];                // # of steps between prints
    string l1_str             = argv[4];                // attraction strength (must be < l2)
    string l2_str             = argv[5];                // attraction range (must be > l1)
    string v0_str             = argv[6];                // tumor cell crawling speed
    string Dr0_str             = argv[7];                // initial angular diffusion
    string aCalA0_str             = argv[8];                // spread of velocity around cell perimeter
    string kecm_str         = argv[9];                // ecm adhesion strength
    string ecmbreak_str     = argv[10];                // ecm adhesion range
    string kv_str           = argv[11];             // ka
    string ka_str           = argv[12];             // kl
    string kc_str           = argv[13];             // kc
    string kb_str           = argv[14];             // kb
    string M_str            = argv[15];
    string P0_str           = argv[16];
    string volumeRatio_str  = argv[17];
    string seed_str         = argv[18];                // seed for rng
    string positionFile     = argv[19];                // output file string
    
    /*
    // read in parameters from command line input
    string inputFile        = "/Users/yitongzheng/Documents/Corey/tumor3D/22div.test";                // input file with initial configuration
    string NT_str           = "100000";                // # of time steps
    string NPRINTSKIP_str   = "1000";                // # of steps between prints
    string l1_str           = "0.0";                // attraction strength (must be < l2)
    string l2_str           = "0.0";                // attraction range (must be > l1)
    string v0_str           = "0.0001";                // tumor cell crawling speed
    string Dr0_str          = "0.01";                // initial angular diffusion
    string aCalA0_str       = "1.1";                // calA0 of adipocyte
    string kecm_str         = "0.0";                // ecm adhesion strength
    string ecmbreak_str     = "0.0";                // ecm adhesion range
    string kv_str           = "10.0";             // kv
    string ka_str           = "0.1";             // ka
    string kc_str           = "0.03";             // kc
    string kb_str           = "0.0";             // kb
    string M_str            = "1.0";
    string P0_str           = "0.0001";
    string volumeRatio_str  = "45";                //
    string seed_str         = "17";                // seed for rng
    string positionFile     = "/Users/yitongzheng/Documents/Corey/tumor3D/P.pos";                // output file string
     */

    // using sstreams to get parameters
    stringstream NTss(NT_str);
    stringstream NPRINTSKIPss(NPRINTSKIP_str);
    stringstream l1ss(l1_str);
    stringstream l2ss(l2_str);
    stringstream v0ss(v0_str);
    stringstream Dr0ss(Dr0_str);
    stringstream aCalA0ss(aCalA0_str);
    stringstream kecmss(kecm_str);
    stringstream ecmbreakss(ecmbreak_str);
    stringstream kvss(kv_str);
    stringstream kass(ka_str);
    stringstream kcss(kc_str);
    stringstream kbss(kb_str);
    stringstream P0ss(P0_str);
    stringstream Mss(M_str);
    stringstream volumeRatioss(volumeRatio_str);
    stringstream seedss(seed_str);

    // read into data
    NTss             >> NTdbl;
    NPRINTSKIPss     >> NPRINTSKIPdbl;
    l1ss             >> l1;
    l2ss             >> l2;
    v0ss             >> v0;
    Dr0ss             >> Dr0;
    aCalA0ss             >> aCalA0;
    kecmss            >> kecm;
    ecmbreakss         >> ecmbreak;
    kvss             >> kv;
    kass             >> ka;
    kcss            >> kc;
    kbss            >> kb;
    Mss             >> M;
    P0ss            >> P0;
    volumeRatioss     >> volumeRatio;
    seedss             >> seed;

    // cast step dbls to ints
    NT = (int)NTdbl;
    NPRINTSKIP = (int)NPRINTSKIPdbl;

    // instantiate object
    tumor3D tumor3Dobj(inputFile,seed);

    // open position config file
    tumor3Dobj.openPosObject(positionFile);

    // set spring constants
    tumor3Dobj.setkv(kv);
    tumor3Dobj.setka(ka);
    tumor3Dobj.setkb(kb);
    tumor3Dobj.setkc(kc);
    tumor3Dobj.setgamtt(gamtt);

    // activity parameters
    tumor3Dobj.setv0(v0);
    tumor3Dobj.setDr0(Dr0);
    tumor3Dobj.setDs(0.1);
    tumor3Dobj.setkecm(kecm);
    tumor3Dobj.setecmbreak(ecmbreak);

    // attraction parameters
    tumor3Dobj.setl1(l1);
    tumor3Dobj.setl2(l2);

    // time step in MD time units
    tumor3Dobj.setdt(dt0);

    //read in standard polyhedron 42 vertices and adjust based on aCalA0
    tumor3Dobj.readPolyhedron();
    tumor3Dobj.initializePolyhedron(aCalA0);
    
    // initialize neighbor linked list
    tumor3Dobj.initializeNeighborLinkedList3D(boxLengthScale);

    // invasion
    cout.precision(10);
    cout << "Running invasion protocol..." << endl;
    // tumor3Dobj.invasion(invasionForceUpdate,dDr,dPsi,Drmin,NT,NPRINTSKIP);
    dDr = 0.01;
    dPsi = 0.01;
    Drmin = 0.01;

    tumor3Dobj.invasionConstP(invasionForceUpdate,M,P0,0,dDr,dPsi,Drmin,NT,NPRINTSKIP);

    // say goodbye
    cout << "\n** Finished interfaceInvasion.cpp, ending. " << endl;

    return 0;
}
 


