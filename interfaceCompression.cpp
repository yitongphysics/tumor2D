// header files
/*
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
	string volumeRatio_str 	= "64.0";
    string aspectRatio_str  = "1.0";
	string prt_str 			= "0.2";
    string kv_str           = "0.01";             // ka
    string ka_str           = "0.01";             // kl
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
	// using sstreams to get p
/*arameters
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
 
 
 
 
 
 // initialize neighbor linked list
 void dpm::initializeNeighborLinkedList3D(double boxLengthScale) {
     // local variables
     double llscale;
     int i, d, scx, scy, scz, boxid;
     int ci, gi, vi;
     
     // print to console
     //cout << "** initializing neighbor linked list, boxLengthScale = " << boxLengthScale;

     // get largest radius as llscale
     llscale = 2.0 * (*max_element(r.begin(),r.end()));

     // initialize box length vectors
     NBX = 1;
     sb.resize(NDIM);
     lb.resize(NDIM);
     for (d = 0; d < NDIM; d++) {
         // determine number of cells along given dimension by rmax
         sb[d] = floor(L[d] / (boxLengthScale * llscale));

         // just in case, if < 3, change to 3 so box neighbor checking will work
         if (sb[d] < 3)
             sb[d] = 3;

         // determine box length by number of cells
         lb[d] = L[d] / sb[d];

         // count total number of cells
         NBX *= sb[d];
     }

     // initialize list of box nearest neighbors
     scx = sb[0];
     scy = sb[1];
     scz = sb[2];
     NBX = scx * scy * scz;
     NNN = 13; // number of nearest neighbors for each cell in 3D

     nn.resize(NBX);

     // loop over cells, save forward neighbors for each box
     for (i = 0; i < NBX; i++) {
         // reshape entry
         nn[i].resize(NNN);

         int ix = i % scx;
         int iy = (i / scx) % scy;
         int iz = i / (scx * scy);

         // Index for nn[i]
         int nn_idx = 0;

         for (int dz = 0; dz <= 1; ++dz) {
             for (int dy = -1; dy <= 1; ++dy) {
                 for (int dx = -1; dx <= 1; ++dx) {
                     // Skip the cases that lead to double counting
                     if (dz == 0 && dy == 0 && dx <= 0) continue;
                     if (dz == 0 && dy == -1) continue;
                     if (dz == 0 && dy == 0 && dx == -1) continue;
                     
                     int nx = ix + dx;
                     int ny = iy + dy;
                     int nz = iz + dz;

                     // Apply periodic boundary conditions if enabled
                     if (pbc[0]) nx = (nx + scx) % scx;
                     if (pbc[1]) ny = (ny + scy) % scy;
                     if (pbc[2]) nz = (nz + scz) % scz;

                     // Check if the neighbor cell is within the simulation box (non-PBC case)
                     bool inside_box = (nx >= 0 && nx < scx) && (ny >= 0 && ny < scy) && (nz >= 0 && nz < scz);

                     // Calculate linear neighbor index and store it in the nn vector
                     int neighbor_idx = inside_box ? (nx + ny * scx + nz * scx * scy) : -1;
                     nn[i][nn_idx] = neighbor_idx;

                     ++nn_idx;
                 }
             }
         }
     }


     // linked-list variables
     head.resize(NBX);
     last.resize(NBX);
     list.resize(NVTOT + 1);
     C_list.resize(NVTOT);

     // print box info to console
     //cout << ";  initially NBX = " << NBX << " ..." << endl;

     fill(list.begin(), list.end(), 0);
     fill(head.begin(), head.end(), 0);
     fill(last.begin(), last.end(), 0);
     fill(C_list.begin(), C_list.end(), 0);
     
     for (gi = 0; gi < NVTOT; gi++) {
         int ix = int(x[NDIM * gi] / L[0] * scx);
         int iy = int(x[NDIM * gi+1] / L[1] * scy);
         int iz = int(x[NDIM * gi+2] / L[2] * scz);

         // Ensure the particle is within the simulation box bounds
         if (ix >= scx) ix = scx - 1;
         if (iy >= scy) iy = scy - 1;
         if (iz >= scz) iz = scz - 1;
         if (ix < 0) ix = 0;
         if (iy < 0) iy = 0;
         if (iz < 0) iz = 0;

         boxid = ix + iy * scx + iz * scx * scy;
         if (head[boxid] == 0) {
             head[boxid] = gi + 1;
             last[boxid] = gi + 1;
         }
         else {
             list[last[boxid]] = gi + 1;
             last[boxid] = gi + 1;
         }
         
         //initialize C_list
         cindices(ci,vi, gi);
         C_list[gi] = ci;
     }
     /*
     NBX=1;
     fill(list.begin(), list.end(), 0);
     fill(head.begin(), head.end(), 0);
     fill(last.begin(), last.end(), 0);
     for (int gi = 0; gi < NVTOT; gi++) {
         if (head[0] == 0) {
             head[0] = gi + 1;
             last[0] = gi + 1;
         }
         else {
             list[last[0]] = gi + 1;
             last[0] = gi + 1;
         }
     }
     list[0] = 1;
      *//*
 }
*/
