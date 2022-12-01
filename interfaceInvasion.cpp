// File to read in initialized interface configuration, run invasion protocol
// 
// Reads in configuration, sets constants, invades
// 
// Compilation command:
// g++ -O3 --std=c++11 -I src main/tumor2D/interfaceInvasion.cpp src/*.cpp -o tumor.o
// 
// Example execution:
// ./tumor.o tumor_input.test 1e5 1e3 0.01 0.02 0.1 0.1 0.2 0.5 0.1 0.9 0.01 1e-3 1 pos.test
// 
// 


// header files
#include "tumor2D.h"
#include <sstream>

// preprocessor macros
#define NDIM 2

// namspace
using namespace std;

// global constants
//const double ka = 1.0;					// area force spring constant (should be unit)
//const double kl = 0.001; 				// contractility spring constant
//const double kc = 0.2;                    // interaction force spring constant
//const double kb = 0.1;					// bending energy
const double gamtt = 0.0; 				// surface tension
const double boxLengthScale = 2;		// neighbor list box size in units of initial l0
const double dt0 = 10;				// initial magnitude of time step in units of MD time
const double Ftol = 1e-5;

int main(int argc, char const *argv[])
{
	// forces for simulations
	tumor2DMemFn invasionForceUpdate = nullptr;
	invasionForceUpdate = &tumor2D::stickyTumorInterfaceForceUpdate;

	// local variables to be read in
	int NT, NPRINTSKIP, seed;
	double NTdbl, NPRINTSKIPdbl, l1, l2, v0, Dr0, Ds, kecm, ecmbreak, dDr, dPsi, Drmin, ka, kl, kb, kc, M, P0, g0;

	// read in parameters from command line input
	string inputFile 		= argv[1];				// input file with initial configuration
	string NT_str 			= argv[2];				// # of time steps
	string NPRINTSKIP_str 	= argv[3];				// # of steps between prints
	string l1_str 			= argv[4];				// attraction strength (must be < l2)
	string l2_str 			= argv[5];				// attraction range (must be > l1)
	string v0_str 			= argv[6];				// tumor cell crawling speed
	string Dr0_str 			= argv[7];				// initial angular diffusion
	string Ds_str 			= argv[8];				// spread of velocity around cell perimeter
	string kecm_str 		= argv[9];				// ecm adhesion strength
	string ecmbreak_str 	= argv[10];				// ecm adhesion range
    string ka_str           = argv[11];             // ka
    string kl_str           = argv[12];             // kl
    string kc_str           = argv[13];             // kc
    string kb_str           = argv[14];             // kb
    string M_str            = argv[15];
    string P0_str           = argv[16];
	string g0_str   		= argv[17];				//
	string seed_str 		= argv[18];				// seed for rng
	string positionFile 	= argv[19];				// output file string

	// using sstreams to get parameters
	stringstream NTss(NT_str);
	stringstream NPRINTSKIPss(NPRINTSKIP_str);
	stringstream l1ss(l1_str);
	stringstream l2ss(l2_str);
	stringstream v0ss(v0_str);
	stringstream Dr0ss(Dr0_str);
	stringstream Dsss(Ds_str);
	stringstream kecmss(kecm_str);
	stringstream ecmbreakss(ecmbreak_str);
	stringstream kass(ka_str);
	stringstream klss(kl_str);
    stringstream kcss(kc_str);
    stringstream kbss(kb_str);
    stringstream P0ss(P0_str);
    stringstream Mss(M_str);
	stringstream g0ss(g0_str);
	stringstream seedss(seed_str);

	// read into data
	NTss 			>> NTdbl;
	NPRINTSKIPss 	>> NPRINTSKIPdbl;
	l1ss 			>> l1;
	l2ss 			>> l2;
	v0ss 			>> v0;
	Dr0ss 			>> Dr0;
	Dsss 			>> Ds;
	kecmss			>> kecm;
	ecmbreakss 		>> ecmbreak;
	kass 			>> ka;
	klss 			>> kl;
    kcss            >> kc;
    kbss            >> kb;
    Mss             >> M;
    P0ss            >> P0;
	g0ss    		>> g0;
	seedss 			>> seed;

	// cast step dbls to ints
	NT = (int)NTdbl;
	NPRINTSKIP = (int)NPRINTSKIPdbl;

	// instantiate object
	tumor2D tumor2Dobj(inputFile,seed);

	// open position config file
	tumor2Dobj.openPosObject(positionFile);

	// set spring constants
	tumor2Dobj.setka(ka);
	tumor2Dobj.setkl(kl);
	tumor2Dobj.setkb(kb);
	tumor2Dobj.setkc(kc);
	tumor2Dobj.setgamtt(gamtt);

	// activity parameters
	tumor2Dobj.setv0(v0);
	tumor2Dobj.setDr0(Dr0);
	tumor2Dobj.setDs(Ds);
	tumor2Dobj.setkecm(kecm);
	tumor2Dobj.setecmbreak(ecmbreak);

	// attraction parameters
	tumor2Dobj.setl1(l1);
	tumor2Dobj.setl2(l2);

	// time step in MD time units
	tumor2Dobj.setdt(dt0);

	// initialize neighbor linked list
	tumor2Dobj.initializeNeighborLinkedList2D(boxLengthScale);

	// run FIRE to relax forces fully
	//tumor2Dobj.tumorFIRE(invasionForceUpdate,Ftol,0.2*dt0);
	// invasion
    cout.precision(10);
	cout << "Running invasion protocol..." << endl;
	// tumor2Dobj.invasion(invasionForceUpdate,dDr,dPsi,Drmin,NT,NPRINTSKIP);
    dDr = 0.01;
    dPsi = 0.01;
    Drmin = 0.01;

	tumor2Dobj.invasionConstP(invasionForceUpdate,M,P0,g0,dDr,dPsi,Drmin,NT,NPRINTSKIP);

	// say goodbye
	cout << "\n** Finished interfaceInvasion.cpp, ending. " << endl;

	return 0;
}
