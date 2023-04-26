#ifndef TUMOR3D_H
#define TUMOR3D_H

#include "dpm.h"
#include <random>

class tumor3D;
typedef void (tumor3D::*tumor3DMemFn)(void);

class tumor3D : public dpm{
protected:

    // temporary number of tumor cells
    int tN;

    //energy terms
    double Ua,Ul,Ub,Utest;
    // surface tension info
    double gamtt;
    
    // wall pressures
    std::vector<double> wpress;
    double wpos;
    
    // motility parameters
    double v0, Dr0, Ds, tau, rho0, rhot;
    std::vector<double> psi;
    std::vector<double> Dr;

    //repulsion between particles
    std::vector<double> F_ij;
    std::vector<int> contactTime;
    // parameters to pin adipocytes
    double kecm, ecmbreak;
    std::vector<double> pinpos;
    std::vector<double> pinattach;
    std::vector<double> ifbroken;

public:

    // Constructors and Destructors
    tumor3D(std::string &inputFileString, int seed);
    tumor3D(int n, int seed) : dpm(n,seed) { tN=n; gamtt=0.0; v0=0.0; Dr0=0.0; Ds=0.0; tau=0.0; rho0=0.0; rhot=0.0; kecm=0.0; ecmbreak=0.0; pbc[0]=0; pbc[1]=0; wpress.resize(2); };                     // tumor-only constructor
    tumor3D(int n, int tNval, int seed) : dpm(n,seed) { tN=tNval; gamtt=0.0; v0=0.0; Dr0=0.0; Ds=0.0; tau=0.0; rho0=0.0; rhot=0.0; kecm=0.0; ecmbreak=0.0; pbc[0]=0; pbc[1]=1; wpress.resize(2); };     // tumor+adipocyte constructor (x fixed, y periodic)

    void readPolyhedron();
    void initializePolyhedron(double aCalA0);
    
    // setters
    void setgamtt(double val) { gamtt = val; };
    void setv0(double val) { v0 = val; };
    void setDr0(double val) { Dr0 = val; };
    void setDs(double val) { Ds = val; };
    void settau(double val) { tau = val; };
    void setkecm(double val) { kecm = val; };
    void setecmbreak(double val) { ecmbreak = val; };
    void reCellList3D(double boxLengthScale);
    void reNeighborList3D();

    // initialization
    void initializeTumorInterface(double aCalA0, double volumeRatio, int aNV, int tNV);
    void initializeTumorInterfacePositions(double phi0, double Ftol, double prt, double aspectRatio);

    // biology functions
    void psiDiffusion();
    void crawlerUpdate();

    // force updates
    void resetForcesAndEnergy();
    void tumorShapeForces();
    void repulsiveTumorInterfaceForces();
    void repulsiveTumorInterfaceForceUpdate();
    void stickyTumorInterfaceForces();
    void stickyTumorInterfaceForceUpdate();

    // integrators
    void tumorFIRE(tumor3DMemFn forceCall, double Ftol, double dt0, int div);

    // protocols
    void setupCheck();
    void tumorCompression(double Ftol, double Ptol, double dt0, double dphi0);
    void invasionConstP(tumor3DMemFn forceCall, double M, double P0, double g0, double dDr, double dPsi, double Drmin, int NT, int NPRINTSKIP);

    // print functions
    void printTumorInterface(double t);
};

#endif

