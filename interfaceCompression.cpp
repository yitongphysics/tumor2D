// header files

#include "tumor3D.h"
#include <sstream>
#include <numeric>
// preprocessor macros
#define NDIM 3

// namspace
using namespace std;

// global constants
const double dphi0= 0.01;	   		// packing fraction increment during initial growth step
const double boxLengthScale = 1.5; 	// neighbor list box size
const double phi0 = 0.01;		   	// initial packing fraction
const double dt0 = 0.03;		   	// initial magnitude of time step in units of MD time
const double Ftol = 1e-7;			// force tolerance during energy min

int main(int argc, char const *argv[])
{
	// local variables to be read in
	int NCELLS, aN, tN, aNV, tNV, seed;
	double aDisp, tDisp, aCalA0, tCalA0, volumeRatio, prt, kv, ka, kb, kc, P0, aspectRatio;

	// read in parameters from command line input
 
 
	string aN_str 			= "28";
	string aNV_str 			= "42";
	string tNV_str 			= "1";
	string aDisp_str 		= "0.0";
	string tDisp_str 		= "0.0";
	string aCalA0_str 		= "1.1";//1.02405832074679303822506426513428
	string tCalA0_str 		= "1.0";
	string volumeRatio_str 	= "45.0";
    string aspectRatio_str  = "1.0";
	string prt_str 			= "0.2";
    string kv_str           = "10.0";             // ka
    string ka_str           = "0.1";             // kl
    string kc_str           = "0.03";             // kc
    string kb_str           = "0.0";             // kb
    string P0_str           = "0.0001";
	string seed_str 		= "1";
	string positionFile 	= "/Users/yitongzheng/Documents/Corey/tumor3D/C.test";
     
    /*
    string aN_str             = argv[1];
    string aNV_str             = argv[2];
    string tNV_str             = argv[3];
    string aDisp_str         = argv[4];
    string tDisp_str         = argv[5];
    string aCalA0_str         = argv[6];
    string tCalA0_str         = argv[7];
    string volumeRatio_str     = argv[8];
    string aspectRatio_str  = argv[9];
    string prt_str             = argv[10];
    string kv_str           = argv[11];             // ka
    string ka_str           = argv[12];             // kl
    string kc_str           = argv[13];             // kc
    string kb_str           = argv[14];             // kb
    string P0_str           = argv[15];
    string seed_str         = argv[16];
    string positionFile     = argv[17];
    */
	// using sstreams to get parameters
	stringstream aNss(aN_str);
	stringstream aNVss(aNV_str);
	stringstream tNVss(tNV_str);
	stringstream aDispss(aDisp_str);
	stringstream tDispss(tDisp_str);
	stringstream aCalA0ss(aCalA0_str);
	stringstream tCalA0ss(tCalA0_str);
	stringstream volumeRatioss(volumeRatio_str);
    stringstream aspectRatioss(aspectRatio_str);
	stringstream prtss(prt_str);
    stringstream kvss(kv_str);
    stringstream kass(ka_str);
    stringstream kcss(kc_str);
    stringstream kbss(kb_str);
    stringstream P0ss(P0_str);
	stringstream seedss(seed_str);

	// read into data
	aNss 			>> aN;
	aNVss 			>> aNV;
	tNVss 			>> tNV;
	aDispss 		>> aDisp;
	tDispss 		>> tDisp;
	aCalA0ss 		>> aCalA0;
	tCalA0ss 		>> tCalA0;
	volumeRatioss 	>> volumeRatio;
    aspectRatioss   >> aspectRatio;
	prtss 			>> prt;
    kvss            >> kv;
    kass            >> ka;
    kcss            >> kc;
    kbss            >> kb;
    P0ss            >> P0;
	seedss 			>> seed;
    
	// determine number of tumor cells based on volumeRatio and prt
	//tN = round(aN * volumeRatio * (prt/(1.0 - prt)));
    tN = 1500;
	NCELLS = tN + aN;

    cout.precision(10);
	// instantiate object
	tumor3D tumor3Dobj(NCELLS, tN, seed);
    
	// open position config file
	tumor3Dobj.openPosObject(positionFile);

	// set spring constants
	tumor3Dobj.setkv(kv);
	tumor3Dobj.setka(ka);
	tumor3Dobj.setkb(kb);
	tumor3Dobj.setkc(kc);

    tumor3Dobj.readPolyhedron();
    
	// initialize adipocyte and tumor cells
	tumor3Dobj.initializeTumorInterface(aCalA0, volumeRatio, aNV, tNV);

	// initialize particle positions
	tumor3Dobj.initializeTumorInterfacePositions(phi0, Ftol, prt, aspectRatio);

	// initialize neighbor linked list
	tumor3Dobj.initializeNeighborLinkedList3D(boxLengthScale);

	// -- compression to initial condition
	tumor3Dobj.tumorCompression(Ftol,P0,dt0,dphi0);
    
	// print interface
	tumor3Dobj.printTumorInterface(0.0);


	// end
	cout << "interfaceCompression.cpp completed, interface printed to " << positionFile << ", ending. " << endl;
	return 0;
}

