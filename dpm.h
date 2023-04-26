#ifndef DPM_H
#define DPM_H
#define SNUM         42 // total number of sites/vertices on one cell
#define FNUM         80 // total number of faces on one cell
#define ENUM         120 // total number of edges on one cell
#define NN_NUM         741


#include <iostream>
#include <iomanip>
#include <string>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <functional>
#include <algorithm>
#include <random>


// pointer-to-member function call macro
#define CALL_MEMBER_FN(object, ptrToMember) ((object).*(ptrToMember))

class dpm;
typedef void (dpm::*dpmMemFn)(void);

// global constants
const double PI = 4.0 * atan(1.0);
const int nvmin = 12;

// printing constants
const int w = 10;
const int wnum = 25;
const int pnum = 14;

// FIRE constants
const double alpha0 = 0.15;
const double finc = 1.1;
const double fdec = 0.5;
const double falpha = 0.99;

const int NSKIP = 5000;
const int NMIN = 10;
const int NNEGMAX = 1000;
const int NDELAY = 20;
const int itmax = 5e7;

class dpm
{
protected:
	// int scalars
	int NCELLS;
	int NDIM;
	int NNN;
	int NVTOT;
	int vertDOF;

	// time step size
	double dt;

	// potential energy
	double U;
    double Udpm;
    
	// particle spring constants
	double kv;
	double ka;
	double kb;
	double kc;

	// particle attraction constants
	double l1, l2;

	// boundary parameters
	std::vector<double> L;
	std::vector<bool> pbc;

	// particle shape parameters
	std::vector<double> r;
    double D0_inter_unit;
    
	// indexing variables
	std::vector<int> nv;
    std::vector<double> V0;
	std::vector<int> szList;
    std::vector<int> C_list;
    
    std::vector<std::vector<int>>        edgelist;
    std::vector<std::vector<int>>        f_unit;
    std::vector<std::vector<int>>        nnlist;
    std::vector<double>             A0_unit;
    std::vector<double>             A0;
    std::vector<double>             xyz_unit;

	// dynamical variables
	std::vector<double> x;
	std::vector<double> v;
	std::vector<double> F;

	// Box linked-list variables
	int NBX;
	std::vector<int> sb;
	std::vector<double> lb;
	std::vector<std::vector<int>> nn;
	std::vector<int> head;
	std::vector<int> last;
	std::vector<int> list;

	// output objects
	std::ofstream posout;

public:
	// Constructors and Destructors
	dpm(int ndim);
	dpm(int n, int ndim, int seed);
	dpm(int n, int seed) : dpm(n, 3, seed) {}
	~dpm();

	// -- G E T T E R S
	// cell shape indexing + information
	int gindex(int ci, int vi);
	void cindices(int &ci, int &vi, int gi);
	double volume(int ci);
	double area(int ci);
	void com3D(int ci, double &cx, double &cy, double &cz);
    void cof3D(int ci, double &fx, double &fy, double &fz);
    void cov3D(int ci, double &vx, double &vy, double &vz);
    double NN3D(int gi);
	double vertexPackingFraction3D();
	double vertexPreferredPackingFraction3D();

	// Setters
	void setdt(double val);
    void setka(double val) { ka = val; };
	void setkv(double val) { kv = val; };
	void setkb(double val) { kb = val; };
	void setkc(double val) { kc = val; };
	void setl1(double val) { l1 = val; };
	void setl2(double val) { l2 = val; };

	// File openers
	void openPosObject(std::string &str)
	{
		posout.open(str.c_str());
		if (!posout.is_open()) {
			std::cerr << "	ERROR: posout could not open " << str << "..." << std::endl;
			exit(1);
		}
		else
			std::cout << "** Opening pos file " << str << " ..." << std::endl;
	}
    
	// Initialize particles (three dimensions)
	void initializeVertexShapeParameters(int ci, double calA0, double lenscale);
	void initializeNeighborLinkedList3D(double boxLengthScale);

	// editing & updating
	void scaleParticleSizes3D(double scaleFactor);
};

#endif
