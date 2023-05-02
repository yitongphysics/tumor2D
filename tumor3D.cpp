#include "tumor3D.h"
#include "dpm.h"
#include <numeric>

// namespace
using namespace std;



/*********************************

    C O N S T R U C T O R S

    &

    D E S T R U C T O R S

**********************************/


// read-in constructor
tumor3D::tumor3D(string &inputFileStr,int seed) : dpm(3) {
    // open file
    ifstream inputobj(inputFileStr.c_str());
    if (!inputobj.is_open()){
        cerr << "** ERROR: In tumor3D constructor, could not open file " << inputFileStr << ", ending here. " << endl;
        exit(1);
    }

    // set variables to default
    gamtt=0.0; v0=0.0; Dr0=0.0; Ds=0.0; tau=0.0; rho0=0.0; rhot=0.0; kecm=0.0; ecmbreak=0.0; pbc[0]=0; pbc[1]=1;  pbc[2]=1;

    // local variables
    int nvtmp, ci, vi, i;
    double val;
    double lxtmp, lytmp, wpostmp;
    double s1, s2, s3;
    double wp1, wp2;
    double V0tmp, psitmp, a0tmp;
    double xtmp, ytmp, ztmp, rtmp, vxtmp, vytmp, vztmp;
    string inputStr;

    // LINE 1: should be NEWFR
    getline(inputobj, inputStr);
    if (inputStr.compare(0,5,"NEWFR") != 0){
        cerr << "** ERROR: In tumor3D constructor, first line of input file NOT NEWFR, first line = " << inputStr << ". Ending." << endl;
        exit(1);
    }

    // read in simulation information from header
    getline(inputobj, inputStr);
    sscanf(inputStr.c_str(),"NUMCL %d %d",&NCELLS,&tN);
    //cout << "\t ** " << inputStr << endl;

    // verify input file
    if (NCELLS< 1){
        cerr << "** ERROR: in tumor3D constructor, NCELLStmp = " << NCELLS << ". Ending here." << endl;
        exit(1);
    }
    else if (tN < 1){
        cerr << "** ERROR: in tumor3D constructor, tNtmp = " << tN << ". Ending here." << endl;
        exit(1);
    }

    getline(inputobj, inputStr);
    sscanf(inputStr.c_str(),"TSTEP %lf",&val);
    //cout << "\t ** " << inputStr << endl;

    getline(inputobj, inputStr);
    sscanf(inputStr.c_str(),"TCURR %lf",&val);
    //cout << "\t ** " << inputStr << endl;
    
    getline(inputobj, inputStr);
    sscanf(inputStr.c_str(),"PACKF %lf",&val);
    //cout << "\t ** " << inputStr << endl;

    // initialize box lengths
    getline(inputobj, inputStr);
    sscanf(inputStr.c_str(),"BOXSZ %lf %lf %lf",&lxtmp,&lytmp,&wpostmp);
    //cout << "\t ** " << inputStr << endl;
    /*
    if (wpostmp!=0) {
        if (abs(log10(abs(wpostmp)))>3 || abs(wpostmp)>4) {
            wpostmp = 0;
        }
    }
    */
    L.at(0) = lxtmp;
    L.at(1) = lytmp;
    L.at(2) = lytmp;
    wpos = wpostmp;
    // initialize stress
    getline(inputobj, inputStr);
    sscanf(inputStr.c_str(),"STRSS %lf %lf %lf",&s1,&s2,&s3);
    //cout << "\t ** " << inputStr << endl;

    // initialize wall pressure
    getline(inputobj, inputStr);
    sscanf(inputStr.c_str(),"WPRSS %lf %lf",&wp1,&wp2);
    //cout << "\t ** " << inputStr << endl;

    wpress.resize(2);
    wpress.at(0) = wp1;
    wpress.at(1) = wp2;
    
    // szList and nv (keep track of global vertex indices)
    nv.resize(NCELLS);
    V0.resize(NCELLS);
    szList.resize(NCELLS);

    // initialize ecm + crawling variables
    psi.resize(tN*NDIM);
    Dr.resize(tN);
    
    fill(psi.begin(), psi.end(), 0.0);
    fill(Dr.begin(), Dr.end(), 0.0);
    
    pinpos.resize(NDIM * (NCELLS - tN));
    pinattach.resize(NCELLS - tN);
    ifbroken.resize(NCELLS - tN);

    // initialize NVTOT to 0
    NVTOT = 0;

    // loop over cells, read in coordinates
    cout << "\n\n** LOOPING OVER DATA, PRINTING INFO..." << endl;
    for (ci=0; ci<NCELLS; ci++){
        // first parse cell info
        getline(inputobj, inputStr);
        if (ci < tN){
            sscanf(inputStr.c_str(),"CINFO %d %*d %*d %lf %*lf %lf %*f",&nvtmp,&V0tmp,&psitmp);
            //psi.at(ci) = psitmp;
        }
        else
            sscanf(inputStr.c_str(),"CINFO %d %*d %*d %lf %*lf %*lf %*lf %*lf %*lf",&nvtmp,&V0tmp);

        // print to console
        //cout << "\t ** " << inputStr << endl;

        // store in vectors
        nv.at(ci) = nvtmp;
        V0.at(ci) = V0tmp;
        // update NVTOT
        NVTOT += nvtmp;

        // increment szList
        if (ci > 0)
            szList.at(ci) = szList.at(ci-1) + nv.at(ci-1);

        // loop over vertices, store coordinates
        for (vi=0; vi<nvtmp; vi++){
            // parse vertex coordinate info
            getline(inputobj, inputStr);
            sscanf(inputStr.c_str(),"VINFO %*d %*d %lf %lf %lf %lf %lf %lf %lf %lf",&xtmp,&ytmp,&ztmp,&rtmp,&a0tmp,&vxtmp,&vytmp,&vztmp);

            // push back
            x.push_back(xtmp);
            x.push_back(ytmp);
            x.push_back(ztmp);
            r.push_back(rtmp);
            v.push_back(vxtmp);
            v.push_back(vytmp);
            v.push_back(vztmp);
        }
    }
    vertDOF = NDIM * NVTOT;
    x.resize(vertDOF);
    v.resize(vertDOF);
    F.resize(vertDOF);
    
    // LINE last: should be ENDFR
    getline(inputobj, inputStr);
    if (inputStr.compare(0,5,"ENDFR") != 0){
        cerr << "** ERROR: In tumor3D constructor, last line of input file NOT ENDFR, first line = " << inputStr << ". Ending." << endl;
        exit(1);
    }
    cout << "** FINISHED LOOPING OVER DATA, CLOSING INPUT FILE OBJ\n" << endl;

    //set rho0
    for(ci=tN;ci<NCELLS;ci++){
        rho0 +=V0[ci];
    }
    rho0 = rho0/(NCELLS - tN);
    //set rhot
    rhot=1;
    for(ci=0;ci<tN;ci++){
        if (rhot>V0[ci]) {
            rhot=V0[ci];
        }
    }

    // close input file object
    inputobj.close();

    // seed random number generator
    srand48(seed);
}

void tumor3D::readPolyhedron(){
    // Read in D0, V0, and A0 for the unit cell, and edge list and face list
    int i;
    int hm = 1;
    A0_unit.resize(FNUM, 0.0);
    A0.resize(FNUM, 0.0);
    
    f_unit.resize(FNUM);
    for(i=0; i<FNUM; i++){
        f_unit[i].resize(3,0);
    }
    
    edgelist.resize(ENUM);
    for(i=0; i<ENUM; i++){
        edgelist[i].resize(6,0);
    }
    
    nnlist.resize(NN_NUM);
    for(i=0; i<NN_NUM; i++){
        nnlist[i].resize(2,0);
    }
    
    xyz_unit.resize(3 * SNUM, 0.0);
    
    if(hm==0){
        std::ifstream pfile1("/gpfs/gibbs/project/ohern/yz974/tumor3D/src/DVA.txt");
        if(pfile1.good())
        {
            pfile1 >> D0_inter_unit;
            for(int i = 0; i < FNUM; i++)
                pfile1 >> A0_unit[i];
        }
        else{
            pfile1.close();
            cout << "Could not read DVA.txt" <<endl;
            exit(-1);
        }
        pfile1.close();
        
        std::ifstream pfile2("/gpfs/gibbs/project/ohern/yz974/tumor3D/src/F_unit.txt");
        if(pfile2.good())
        {
            for (int i = 0; i < FNUM; i++)
                for (int j = 0; j < 3; j++)
                    pfile2 >> f_unit[i][j];
        }
        else{
            pfile2.close();
            cout << "Could not read F_unit.txt" <<endl;
            exit(-1);
        }
        pfile2.close();
        
        std::ifstream pfile3("/gpfs/gibbs/project/ohern/yz974/tumor3D/src/EdgeList.txt");
        if(pfile3.good())
        {
            for (int i = 0; i < ENUM; i++)
                for (int j = 0; j < 6; j++)
                    pfile3 >> edgelist[i][j];
        }
        else{
            pfile3.close();
            cout << "Could not read EdgeList.txt" <<endl;
            exit(-1);
        }
        pfile3.close();
        
        std::ifstream pfile4("/gpfs/gibbs/project/ohern/yz974/tumor3D/src/NonNeighborList.txt");
        if(pfile4.good())
        {
            for (int i = 0; i < NN_NUM; i++)
                for (int j = 0; j < 2; j++)
                    pfile4 >> nnlist[i][j];
        }
        else{
            pfile4.close();
            cout << "Could not read NN.txt" <<endl;
            exit(-1);
        }
        pfile4.close();
        
        std::ifstream pfile5("/gpfs/gibbs/project/ohern/yz974/tumor3D/src/XYZ_unit.txt");
        if(pfile5.good())
        {
            for(int i = 0; i < 3 * SNUM; i++)
                pfile5 >> xyz_unit[i];
        }
        else{
            pfile5.close();
            cout << "Could not read XYZ.txt" <<endl;
            exit(-1);
        }
        pfile5.close();
    }
    else{
        std::ifstream pfile1("/Users/yitongzheng/Documents/Corey/tumor3D/src/DVA.txt");
        if(pfile1.good())
        {
            pfile1 >> D0_inter_unit;
            for(int i = 0; i < FNUM; i++)
                pfile1 >> A0_unit[i];
        }
        else{
            pfile1.close();
            cout << "Could not read DVA.txt" <<endl;
            exit(-1);
        }
        pfile1.close();
        
        std::ifstream pfile2("/Users/yitongzheng/Documents/Corey/tumor3D/src/F_unit.txt");
        if(pfile2.good())
        {
            for (int i = 0; i < FNUM; i++)
                for (int j = 0; j < 3; j++)
                    pfile2 >> f_unit[i][j];
        }
        else{
            pfile2.close();
            cout << "Could not read F_unit.txt" <<endl;
            exit(-1);
        }
        pfile2.close();
        
        std::ifstream pfile3("/Users/yitongzheng/Documents/Corey/tumor3D/src/EdgeList.txt");
        if(pfile3.good())
        {
            for (int i = 0; i < ENUM; i++)
                for (int j = 0; j < 6; j++)
                    pfile3 >> edgelist[i][j];
        }
        else{
            pfile3.close();
            cout << "Could not read EdgeList.txt" <<endl;
            exit(-1);
        }
        pfile3.close();
        
        std::ifstream pfile4("/Users/yitongzheng/Documents/Corey/tumor3D/src/NonNeighborList.txt");
        if(pfile4.good())
        {
            for (int i = 0; i < NN_NUM; i++)
                for (int j = 0; j < 2; j++)
                    pfile4 >> nnlist[i][j];
        }
        else{
            pfile4.close();
            cout << "Could not read NN.txt" <<endl;
            exit(-1);
        }
        pfile4.close();
        
        std::ifstream pfile5("/Users/yitongzheng/Documents/Corey/tumor3D/src/XYZ_unit.txt");
        if(pfile5.good())
        {
            for(int i = 0; i < 3 * SNUM; i++)
                pfile5 >> xyz_unit[i];
        }
        else{
            pfile5.close();
            cout << "Could not read XYZ.txt" <<endl;
            exit(-1);
        }
        pfile5.close();
    }
    
}


/*********************************

    T U M O R   C E L L

    I N I T I A L I Z A T I O N

**********************************/
// initialize neighbor linked list: with wall position change
void tumor3D::reCellList3D(double boxLengthScale) {
    // local variables
    double llscale, Lx = L[0] - wpos;
    int i, d, scx, scy, scz;
    // print to console
    //cout << "** initializing neighbor linked list, boxLengthScale = " << boxLengthScale;
    
    // get largest radius as llscale
    llscale = 2.0 * (*max_element(r.begin(),r.end()));
    
    // initialize box length vectors
    NBX = 1;
    for (d = 0; d < NDIM; d++) {
        if(d==0){
            // determine number of cells along given dimension by rmax
            sb[d] = floor(Lx/ (boxLengthScale * llscale));
            
            // just in case, if < 3, change to 3 so box neighbor checking will work
            if (sb[d] < 3)
                sb[d] = 3;
            
            // determine box length by number of cells
            lb[d] = Lx/ sb[d];
        }
        else{
            // determine number of cells along given dimension by rmax
            sb[d] = floor(L[d] / (boxLengthScale * llscale));
            
            // just in case, if < 3, change to 3 so box neighbor checking will work
            if (sb[d] < 3)
                sb[d] = 3;
            
            // determine box length by number of cells
            lb[d] = L[d]/ sb[d];
        }
        
        // count total number of cells
        NBX *= sb[d];
    }
    
    // initialize list of box nearest neighbors
    scx = sb[0];
    scy = sb[1];
    scz = sb[2];
    NBX = scx * scy * scz;
    
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
                    
                    int nx = ix + dx;
                    int ny = iy + dy;
                    int nz = iz + dz;
                    
                    // Apply periodic boundary conditions if enabled
                    ny = (ny + scy) % scy;
                    nz = (nz + scz) % scz;
                    
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
}
void tumor3D::reNeighborList3D() {
    double Lx = L[0] - wpos;
    int scx, scy, scz, boxid;
    double L1_inv = 1.0 / L[1];
    double L2_inv = 1.0 / L[2];
    double Lx_inv = 1.0 / Lx;
    
    fill(list.begin(), list.end(), 0);
    fill(head.begin(), head.end(), 0);
    fill(last.begin(), last.end(), 0);
    
    scx = sb[0];
    scy = sb[1];
    scz = sb[2];
    for (int gi = 0; gi < NVTOT; gi++) {
        int ix = int((x[NDIM * gi]-wpos) * Lx_inv * scx);
        int iy = int(x[NDIM * gi+1] * L1_inv * scy);
        int iz = int(x[NDIM * gi+2] * L2_inv * scz);

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
    */
    
}

void tumor3D::initializePolyhedron(double aCalA0){
    int fi;
    double A_sum, Rc;
    double calA = 1.02405832074679303822506426513428;
    
    A_sum=0.0;
    for (fi=0; fi<FNUM; fi++){
        A_sum += A0_unit[fi];
    }
    A_sum = 11.6659313917183045106185090844519;
    
    //V0=1
    Rc = pow((  pow(A_sum, 1.5)/sqrt(PI)/6.0/calA  ), -1.0 / 3.0);
    //adjust A0 based on aCalA0
    Rc *= pow(aCalA0/calA,1.0/3.0);
    for (fi=0; fi<FNUM; fi++){
        A0[fi] = A0_unit[fi] * Rc * Rc;
    }
    
    //if (V0[tN+1] != 1.0) {
      //  cout << "Preferred volume is not 1.0. Exit.";
        //exit(-1);
    //}
    
}
void tumor3D::initializeTumorInterface(double aCalA0, double volumeRatio, int aNV, int tNV){
    // local variables
    int ci, nvtmp, fi;
    double lenscale, Rc, A_sum;
    double calA = 1.02405832074679303822506426513428;
    // print to console
    cout << "** initializing tumor and adipocyte DPM particles in 3D..." << endl;
    cout << "** setting up nv + szList, setting shape parameters and initializing indexing ..." << endl;

    // szList and nv (keep track of global vertex indices)
    nv.resize(NCELLS);
    V0.resize(NCELLS);
    szList.resize(NCELLS);

    // initialize ecm + crawling variables
    psi.resize(tN*NDIM);
    Dr.resize(tN);
    
    fill(psi.begin(), psi.end(), 0.0);
    fill(Dr.begin(), Dr.end(), 0.0);

    pinpos.resize(NDIM * (NCELLS - tN));
    pinattach.resize(NCELLS - tN);
    ifbroken.resize(NCELLS - tN);

    // initialize number of vertices on each cell
    nv.at(0) = tNV;
    NVTOT = tNV;
    for (ci=1; ci<NCELLS; ci++){
        if (ci < tN){
            nvtmp = tNV;//floor(tDisp*tNV*grv + tNV);
            V0.at(ci) = 1.0/volumeRatio;
        }
        else{
            nvtmp = aNV;
            V0.at(ci) = 1.0;
        }
        // store size of cell ci
        nv.at(ci) = nvtmp;
        szList.at(ci) = szList.at(ci-1) + nv.at(ci-1);

        // add to total NV count
        NVTOT += nvtmp;
        
    }
    vertDOF = NDIM * NVTOT;

    // resize shape paramters
    x.resize(vertDOF);
    v.resize(vertDOF);
    F.resize(vertDOF);
    r.resize(NVTOT);
    
    
    A_sum=0.0;
    for (fi=0; fi<FNUM; fi++){
        A_sum += A0_unit[fi];
    }
    A_sum = 11.6659313917183045106185090844519;
    
    //V0=1
    Rc = pow((  pow(A_sum, 1.5)/sqrt(PI)/6.0/calA  ), -1.0 / 3.0);
    //adjust A0 based on aCalA0
    Rc *= pow(aCalA0/calA,1.0/3.0);
    for (fi=0; fi<FNUM; fi++){
        A0[fi] = A0_unit[fi] * Rc * Rc;
    }
    
    // initialize particle sizes based on volumeRatio
    for (ci=0; ci<NCELLS; ci++){
        if (ci < tN){
            lenscale = 1.0/pow(volumeRatio, 1.0/3.0);
            
            //bidisperse
            if (ci < tN/2) {
                lenscale = lenscale * 0.9;
            }
            else{
                lenscale = lenscale * 1.1;
            }
            
            initializeVertexShapeParameters(ci,1.0,lenscale);
        }
        else{
            lenscale = Rc;
            initializeVertexShapeParameters(ci,aCalA0,lenscale);
        }
    }
     
    /*
    A_sum=0.0;
        for (fi=0; fi<FNUM; fi++){
            A_sum += A0_unit[fi];
        }
        Rc = pow(aCalA0 / (pow(A_sum, 1.5) / sqrt(PI) / 6.0), 1.0 / 3.0);
        
        for (fi=0; fi<FNUM; fi++){
            A0[fi] = A0_unit[fi] * Rc * Rc;
        }
        
        // initialize particle sizes based on volumeRatio
        for (ci=0; ci<NCELLS; ci++){
            if (ci < tN){
                lenscale = 1.0/pow(volumeRatio, 1.0/3.0);
                
                //bidisperse
                if (ci < tN/2) {
                    lenscale = lenscale * 0.9;
                }
                else{
                    lenscale = lenscale * 1.1;
                }
                
                initializeVertexShapeParameters(ci,1.0,lenscale);
            }
            else{
                lenscale = Rc;
                initializeVertexShapeParameters(ci,aCalA0,lenscale);
            }
        }
    */
}


void tumor3D::initializeTumorInterfacePositions(double phi0, double Ftol, double prt, double aspectRatio){
    // local variables
    int i, d, ci, cj, vi, gi, cellDOF = NDIM * NCELLS;
    double volumeSum, xtra = 2.0, xi, Ldiv, Ltmp;
    double aspect_ratio = aspectRatio;
    // local disk vectors
    vector<double> drad(NCELLS, 0.0);
    vector<double> dpos(cellDOF, 0.0);
    vector<double> dv(cellDOF, 0.0);
    vector<double> dF(cellDOF, 0.0);

    // print to console
    cout << "** initializing particle positions using 3D SP model and FIRE relaxation ..." << endl;

    // initialize box size based on packing fraction
    //if N=1
    volumeSum = 0.0;
    for (ci = 0; ci < NCELLS; ci++){
        if (nv.at(ci)==1) {
            volumeSum += V0[ci];
        } else {
            volumeSum += V0[ci];
        }
    }

    // set box size
    L.at(2) = pow(volumeSum/(aspect_ratio*phi0), 1.0/3.0);
    L.at(1) = pow(volumeSum/(aspect_ratio*phi0), 1.0/3.0);
    L.at(0) = aspect_ratio*L[1];

    // dividing wall position between adipocytes and
    Ldiv = prt * L[0];

    // initialize tumor cell centers in left-hand partition of the box
    for (ci=0; ci<tN; ci++){
        dpos.at(NDIM*ci)         = Ldiv*drand48();
        dpos.at(NDIM*ci + 1)     = L[1]*drand48();
        dpos.at(NDIM*ci + 2)     = L[2]*drand48();
        //randomly
        //dpos.at(NDIM*ci)     = L[0]*drand48();
        //dpos.at(NDIM*ci + 1)     = L[1]*drand48();
    }

    // initialize WAT cell centers to the right
    for (ci=tN; ci<NCELLS; ci++){
        dpos.at(NDIM*ci)         = (L[0] - Ldiv)*drand48() + Ldiv;
        dpos.at(NDIM*ci + 1)     = L[1]*drand48();
        dpos.at(NDIM*ci + 2)     = L[2]*drand48();
        //dpos.at(NDIM*ci)         = ((ci-tN+1)%6)*L[0]/6.0+L[0]/12.0;
        //dpos.at(NDIM*ci + 1)     = ceil((ci-tN+1)/6.0)*L[1]/6.0-L[1]/12.0;
    }
    
    // set radii of SP disks
    //if N=1
    for (ci = 0; ci < NCELLS; ci++){
        if (nv.at(ci)==1) {
            drad.at(ci) = pow(V0.at(ci)/PI*3.0/4.0,1.0/3.0);
        } else {
            drad.at(ci) = pow(V0.at(ci)/PI*3.0/4.0,1.0/3.0) + r.at(ci);
        }
    }
    
    volumeSum = 0.0;
    for (ci=0; ci<NCELLS; ci++){
        volumeSum += 4.0/3.0*PI*drad[ci]*drad[ci]*drad[ci];
    }

    // FIRE VARIABLES
    double P = 0;
    double fnorm = 0;
    double vnorm = 0;
    double alpha = alpha0;

    double dt0 = 1e-2;
    double dtmax = 100 * dt0;
    double dtmin = 1e-8 * dt0;

    int npPos = 0;
    int npNeg = 0;

    int fireit = 0;
    double fcheck = 10 * Ftol;

    // interaction variables
    double rij, sij, dtmp, ftmp, vftmp;
    double dr[NDIM];

    // initial step size
    dt = dt0;

    // loop until force relaxes
    while ((fcheck > Ftol) && fireit < itmax) {
        // FIRE step 1. Compute P
        P = 0.0;
        for (i = 0; i < cellDOF; i++)
            P += dv[i] * dF[i];

        // FIRE step 2. adjust simulation based on net motion of degrees of freedom
        if (P > 0) {
            // increase positive counter
            npPos++;

            // reset negative counter
            npNeg = 0;

            // alter simulation if enough positive steps have been taken
            if (npPos > NMIN) {
                // change time step
                if (dt * finc < dtmax)
                    dt *= finc;

                // decrease alpha
                alpha *= falpha;
            }
        }
        else {
            // reset positive counter
            npPos = 0;

            // increase negative counter
            npNeg++;

            // check if simulation is stuck
            if (npNeg > NNEGMAX){
                cerr << "    ** ERROR: During initial FIRE minimization, P < 0 for too long, so ending." << endl;
                exit(1);
            }

            // take half step backwards, reset velocities
            for (i = 0; i < cellDOF; i++)
            {
                // take half step backwards
                dpos[i] -= 0.5 * dt * dv[i];

                // reset velocities
                dv[i] = 0.0;
            }

            // decrease time step if past initial delay
            if (fireit > NDELAY)
            {
                // decrease time step
                if (dt * fdec > dtmin)
                    dt *= fdec;

                // reset alpha
                alpha = alpha0;
            }
        }

        // FIRE step 3. First VV update
        for (i = 0; i < cellDOF; i++)
            dv[i] += 0.5 * dt * dF[i];

        // FIRE step 4. adjust velocity magnitude
        fnorm = 0.0;
        vnorm = 0.0;
        for (i = 0; i < cellDOF; i++) {
            fnorm += dF[i] * dF[i];
            vnorm += dv[i] * dv[i];
        }
        fnorm = sqrt(fnorm);
        vnorm = sqrt(vnorm);
        if (fnorm > 0) {
            for (i = 0; i < cellDOF; i++)
                dv[i] = (1 - alpha) * dv[i] + alpha * (vnorm / fnorm) * dF[i];
        }

        // FIRE step 4. Second VV update
        for (i = 0; i < cellDOF; i++) {
            dpos[i] += dt * dv[i];
            dF[i] = 0.0;
        }

        // FIRE step 5. Update forces
        for (ci = 0; ci < NCELLS; ci++) {
            for (cj = ci + 1; cj < NCELLS; cj++) {

                // contact distance
                sij = drad[ci] + drad[cj];

                // true distance
                rij = 0.0;
                for (d = 0; d < NDIM; d++) {
                    // get distance element
                    dtmp = dpos[NDIM * cj + d] - dpos[NDIM * ci + d];
                    if (pbc[d])
                        dtmp -= L[d] * round(dtmp / L[d]);

                    // add to true distance
                    rij += dtmp * dtmp;

                    // save in distance array
                    dr[d] = dtmp;
                }
                rij = sqrt(rij);

                // check distances
                if (rij < sij) {
                    // force magnitude
                    ftmp = kc * (1.0 - (rij / sij)) / sij;

                    // add to vectorial force
                    for (d = 0; d < NDIM; d++)
                    {
                        vftmp = ftmp * (dr[d] / rij);
                        dF[NDIM * ci + d] -= vftmp;
                        dF[NDIM * cj + d] += vftmp;
                    }
                }
            }

            // x boundary forces
            xi = dpos[NDIM * ci];
            if (ci < tN) {
                if (xi < drad[ci])
                    dF[NDIM*ci] += (1.0 - (xi/drad[ci]))/drad[ci];
                //warning if (xi > L[0]- drad[ci])
                else if (xi > Ldiv - drad[ci])
                    dF[NDIM*ci] -= (1.0 - ((Ldiv - xi)/drad[ci]))/drad[ci];
                    //dF[NDIM*ci] -= (1.0 - ((L[0] - xi)/drad[ci]))/drad[ci];
            }
            else {
                //warning:if (xi < drad[ci])
                if (xi < Ldiv + drad[ci])
                    dF[NDIM*ci] += (1.0 - ((xi - Ldiv)/drad[ci]))/drad[ci];
                    //dF[NDIM*ci] += (1.0 - (xi/drad[ci]))/drad[ci];
                else if (xi > L[0] - drad[ci])
                    dF[NDIM*ci] -= (1.0 - ((L[0] - xi)/drad[ci]))/drad[ci];
            }
        }

        // FIRE step 5. Final VV update
        for (i = 0; i < cellDOF; i++)
            dF[i] -= 0.000001*dpos[i]*dpos[i];
        
        // FIRE step 5. Final VV update
        for (i = 0; i < cellDOF; i++)
            dv[i] += 0.5 * dt * dF[i];

        // update forces to check
        fcheck = 0.0;
        for (i = 0; i < cellDOF; i++)
            fcheck += dF[i] * dF[i];
        fcheck = sqrt(fcheck / NCELLS);

        //shrink box
        Ltmp =  ((*max_element(dpos.begin(),dpos.end())) + (*max_element(drad.begin(),drad.end())));
        L[0] = L[0] < Ltmp ? L[0] : Ltmp;
        L[1] = L[0];
        L[2] = L[0];
        
        if (volumeSum/L[0]/L[1]/L[2] > 0.2) {
            break;
        }
        // print to console
        if (fireit % NSKIP == 0) {
            cout << endl
                 << endl;
            cout << "===========================================" << endl;
            cout << "        I N I T I A L  S P             " << endl;
            cout << "     F I R E                         " << endl;
            cout << "        M I N I M I Z A T I O N     " << endl;
            cout << "===========================================" << endl;
            cout << endl;
            cout << "    ** fireit = " << fireit << endl;
            cout << "    ** fcheck = " << fcheck << endl;
            cout << "    ** fnorm = " << fnorm << endl;
            cout << "    ** vnorm = " << vnorm << endl;
            cout << "    ** dt = " << dt << endl;
            cout << "    ** P = " << P << endl;
            cout << "    ** Pdir = " << P / (fnorm * vnorm) << endl;
            cout << "    ** alpha = " << alpha << endl;
            cout << "    ** phi0 = " << volumeSum/L[0]/L[1]/L[2] << endl;
            
            for (ci = 0; ci < NCELLS; ci++) {
                for (vi = 0; vi < nv.at(ci); vi++) {
                    // get global vertex index
                    gi = gindex(ci, vi);

                    // length from center to vertex
                    if (nv.at(ci)==1) {
                        x.at(NDIM * gi) = dpos.at(NDIM * ci);
                        x.at(NDIM * gi + 1) = dpos.at(NDIM * ci + 1);
                        x.at(NDIM * gi + 2) = dpos.at(NDIM * ci + 2);
                    } else {
                        
                        dtmp = pow(3.0/4.0/PI*V0.at(ci), 1.0/3.0);
                        
                        // set positions
                        x.at(NDIM * gi) = dtmp*xyz_unit[3*vi] + dpos.at(NDIM * ci);
                        x.at(NDIM * gi + 1) = dtmp*xyz_unit[3*vi+1] + dpos.at(NDIM * ci + 1);
                        x.at(NDIM * gi + 2) = dtmp*xyz_unit[3*vi+2] + dpos.at(NDIM * ci + 2);
                    }
                }
            }
            
            //printTumorInterface(0.0);
        }

        // update iterate
        fireit++;
    }
    // check if FIRE converged
    if (fireit == itmax) {
        cout << "    ** FIRE minimization did not converge, fireit = " << fireit << ", itmax = " << itmax << "; ending." << endl;
        exit(1);
    }
    else {
        cout << endl
             << endl;
        cout << "&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&" << endl;
        cout << "===========================================" << endl;
        cout << "     F I R E                         " << endl;
        cout << "        M I N I M I Z A T I O N     " << endl;
        cout << "    C O N V E R G E D!                 " << endl
             << endl;

        cout << "    (for initial disk minimization) " << endl;
        cout << "===========================================" << endl;
        cout << "&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&" << endl;
        cout << endl;
        cout << "    ** fireit = " << fireit << endl;
        cout << "    ** fcheck = " << fcheck << endl;
        cout << "    ** vnorm = " << vnorm << endl;
        cout << "    ** dt = " << dt << endl;
        cout << "    ** P = " << P << endl;
        cout << "    ** alpha = " << alpha << endl;
    }

    // initialize vertex positions based on cell centers
    //if N=1
    for (ci = 0; ci < NCELLS; ci++) {
        for (vi = 0; vi < nv.at(ci); vi++) {
            // get global vertex index
            gi = gindex(ci, vi);

            // length from center to vertex
            if (nv.at(ci)==1) {
                x.at(NDIM * gi) = dpos.at(NDIM * ci);
                x.at(NDIM * gi + 1) = dpos.at(NDIM * ci + 1);
                x.at(NDIM * gi + 2) = dpos.at(NDIM * ci + 2);
            } else {
                
                dtmp = pow(3.0/4.0/PI*V0.at(ci), 1.0/3.0);
                
                // set positions
                x.at(NDIM * gi) = dtmp*xyz_unit[3*vi] + dpos.at(NDIM * ci);
                x.at(NDIM * gi + 1) = dtmp*xyz_unit[3*vi+1] + dpos.at(NDIM * ci + 1);
                x.at(NDIM * gi + 2) = dtmp*xyz_unit[3*vi+2] + dpos.at(NDIM * ci + 2);
            }
        }
    }
}

/******************************

    B I O L O G I C A L

    F U N C T I O N S

*******************************/

//initialize crawlling direction
void tumor3D::initializePsi(){
    int ci;
    double r1, r2, r3;
    double psi_R;
    default_random_engine gen;
    normal_distribution<double> distribution(0.0, 1.0);
    
    for (ci=0; ci<tN; ci++){
        // generate random variable
        r1 = distribution(gen);
        r2 = distribution(gen);
        r3 = distribution(gen);
        
        psi_R = sqrt(r1*r1+r2*r2+r3*r3);
        
        psi[NDIM*ci    ] = r1/psi_R;
        psi[NDIM*ci + 1] = r2/psi_R;
        psi[NDIM*ci + 2] = r3/psi_R;
    }
}

// update psi based on Dr only
void tumor3D::psiDiffusion(int seed){
    // local variables
    int ci;
    double r1, r2, r3;
    double psi_x, psi_y, psi_z;
    double dpsi_x, dpsi_y, dpsi_z;
    double R_inv;
    
    default_random_engine gen(seed);
    normal_distribution<double> distribution(0.0, 1.0);
    
    double sg = sqrt(Dr0*2.0)*dt;
    // update director for each cell
    for (ci=0; ci<tN; ci++){
        psi_x = psi[ci*NDIM];
        psi_y = psi[ci*NDIM + 1];
        psi_z = psi[ci*NDIM + 2];
        
        // generate random variable
        r1 = distribution(gen);
        r2 = distribution(gen);
        r3 = distribution(gen);
        
        dpsi_x = sg * (psi_y*r3 - psi_z*r2);
        dpsi_y = sg * (psi_z*r1 - psi_x*r3);
        dpsi_z = sg * (psi_x*r2 - psi_y*r1);
        // update director for cell ci
        psi_x += dpsi_x;
        psi_y += dpsi_y;
        psi_z += dpsi_z;
        
        R_inv = 1.0/sqrt(psi_x*psi_x+psi_y*psi_y+psi_z*psi_z);
        
        psi[ci*NDIM    ] = psi_x*R_inv;
        psi[ci*NDIM + 1] = psi_y*R_inv;
        psi[ci*NDIM + 2] = psi_z*R_inv;
    }
}
//tumor cells swim
void tumor3D::crawlerUpdate(){
    // local variables
    int ci;
    // loop over all cells, vertices
    for (ci=0; ci<tN; ci++){
        F[NDIM*ci]     += f0*psi[NDIM*ci];
        F[NDIM*ci + 1] += f0*psi[NDIM*ci + 1];
        F[NDIM*ci + 1] += f0*psi[NDIM*ci + 2];
    }
}


/******************************

    F O R C E

    U P D A T E S

*******************************/


void tumor3D::resetForcesAndEnergy(){
    fill(F.begin(), F.end(), 0.0);
    fill(wpress.begin(), wpress.end(), 0.0);
    U = 0.0;
    Udpm = 0.0;
    Ua=0.0;
    Ul=0.0;
    Ub=0.0;
    Utest=0.0;
}

void tumor3D::tumorShapeForces(){
    // compute DPM forces. Did not include forces between vertices on the same cell.
    int         ns, nc, nf, ne;
    int         ns_start;
    int         s_idx_start;
    //double      KC_half     = kc / 2.0;
    double      cx, cy, cz;
    double L1_inv = 1.0 / L[1];
    double L2_inv = 1.0 / L[2];
    
    vector<double>     px(NVTOT, 0.0), py(NVTOT, 0.0), pz(NVTOT, 0.0);
    vector<double>    fx(NVTOT, 0.0), fy(NVTOT, 0.0), fz(NVTOT, 0.0);

    for (ns = 0; ns < NVTOT; ns++){
        fx[ns] = 0.0;
        fy[ns] = 0.0;
        fz[ns] = 0.0;
    }

    for (ns = 0; ns < NVTOT; ns++){
        ns_start = 3 * ns;
        px[ns] = x[ns_start];
        py[ns] = x[ns_start + 1];
        pz[ns] = x[ns_start + 2];
    }

    vector<double>    A(FNUM, 0.0); // areas of triangles on one DP
    double             V; // volume of one DP
    double             x1, y1, z1;
    double             x2, y2, z2;
    double             x3, y3, z3;
    double             x4, y4, z4;
    double             dx12, dy12, dz12;
    double             dx13, dy13, dz13;
    double             dx14, dy14, dz14;
    double             dx32, dy32, dz32;
    double             dx34, dy34, dz34;
    double             r13, r13_inv;
    double             r_dot_1, r_dot_2, r_dot_3, r_dot_4;
    double             FA2_x, FA2_y, FA2_z;
    double             FA3_x, FA3_y, FA3_z;
    vector<double>     R_cross_32_x(FNUM, 0.0), R_cross_32_y(FNUM, 0.0), R_cross_32_z(FNUM, 0.0);
    vector<double>    R_cross_13_x(FNUM, 0.0), R_cross_13_y(FNUM, 0.0), R_cross_13_z(FNUM, 0.0);
    vector<double>     R_cross_21_x(FNUM, 0.0), R_cross_21_y(FNUM, 0.0), R_cross_21_z(FNUM, 0.0);
    double             dA, dV;
    double             Ax, Ay, Az;
    vector<double>    A_vec_x(FNUM, 0.0), A_vec_y(FNUM, 0.0), A_vec_z(FNUM, 0.0);
    vector<double>     A_rescale_x(FNUM, 0.0), A_rescale_y(FNUM, 0.0), A_rescale_z(FNUM, 0.0);
    double             A_sq;
    double             A_dot; // dot product of two normalized area vector
    double             theta; // face-face angle
    int             i;
    int             ns1, ns2, ns3, ns4;
    int             nf1, nf2;
    vector<int>        face1(3, 0), face2(3, 0);
    //double          D0_inter = D0_inter_unit;
    //double          dd, dnm;
    //double          ftmp, dFx, dFy, dFz;
    
    vector<int> ns1_arr(FNUM);
    vector<int> ns2_arr(FNUM);
    vector<int> ns3_arr(FNUM);
    
    for (nc = tN; nc < NCELLS; nc++){
        s_idx_start = (nc-tN) * SNUM + tN;
        V = 0.0;
        com3D(nc, cx, cy, cz);
        

        for (nf = 0; nf < FNUM; nf++){
            ns1_arr[nf] = f_unit[nf][0] + s_idx_start;
            ns2_arr[nf] = f_unit[nf][1] + s_idx_start;
            ns3_arr[nf] = f_unit[nf][2] + s_idx_start;
        }

        // Area force
        for (nf = 0; nf < FNUM; nf++){
            ns1 = ns1_arr[nf];
            ns2 = ns2_arr[nf];
            ns3 = ns3_arr[nf];

            x1 = px[ns1]-cx; y1 = py[ns1]-cy; z1 = pz[ns1]-cz;
            x2 = px[ns2]-cx; y2 = py[ns2]-cy; z2 = pz[ns2]-cz;
            x3 = px[ns3]-cx; y3 = py[ns3]-cy; z3 = pz[ns3]-cz;
            
            y1 -= L[1] * round(y1 * L1_inv);z1 -= L[2] * round(z1 * L2_inv);
            y2 -= L[1] * round(y2 * L1_inv);z2 -= L[2] * round(z2 * L2_inv);
            y3 -= L[1] * round(y3 * L1_inv);z3 -= L[2] * round(z3 * L2_inv);
            
            // record Ri cross Rj for later volume force calculation
            R_cross_32_x[nf] = y3 * z2 - z3 * y2;
            R_cross_32_y[nf] = z3 * x2 - x3 * z2;
            R_cross_32_z[nf] = x3 * y2 - y3 * x2;
            R_cross_13_x[nf] = y1 * z3 - z1 * y3;
            R_cross_13_y[nf] = z1 * x3 - x1 * z3;
            R_cross_13_z[nf] = x1 * y3 - y1 * x3;
            R_cross_21_x[nf] = y2 * z1 - z2 * y1;
            R_cross_21_y[nf] = z2 * x1 - x2 * z1;
            R_cross_21_z[nf] = x2 * y1 - y2 * x1;

            dx12 = x2 - x1;
            dy12 = y2 - y1;
            dz12 = z2 - z1;
            dx13 = x3 - x1;
            dy13 = y3 - y1;
            dz13 = z3 - z1;

            Ax = (dy12 * dz13 - dz12 * dy13) * 0.5;
            Ay = (dz12 * dx13 - dx12 * dz13) * 0.5;
            Az = (dx12 * dy13 - dy12 * dx13) * 0.5;
            A[nf] = sqrt(Ax * Ax + Ay * Ay + Az * Az);
            A_vec_x[nf] = Ax;
            A_vec_y[nf] = Ay;
            A_vec_z[nf] = Az;
            dA = A[nf] / A0[nf] - 1.0;
            
            U += ka * dA * dA;
            Udpm += ka * dA * dA;
            
            dA = dA * ka / A[nf] / A0[nf];
            FA2_x = dA * (Ay * dz13 - Az * dy13);
            FA2_y = dA * (Az * dx13 - Ax * dz13);
            FA2_z = dA * (Ax * dy13 - Ay * dx13);
            FA3_x = dA * (dy12 * Az - dz12 * Ay);
            FA3_y = dA * (dz12 * Ax - dx12 * Az);
            FA3_z = dA * (dx12 * Ay - dy12 * Ax);

            fx[ns1] -= (FA2_x + FA3_x);
            fy[ns1] -= (FA2_y + FA3_y);
            fz[ns1] -= (FA2_z + FA3_z);
            fx[ns2] += FA2_x;
            fy[ns2] += FA2_y;
            fz[ns2] += FA2_z;
            fx[ns3] += FA3_x;
            fy[ns3] += FA3_y;
            fz[ns3] += FA3_z;

            V += (Ax * x1 + Ay * y1 + Az * z1);
        }
        V /= 3.0;

        // Volume force
        dV = V / V0[nc] - 1.0;
        U += kv * dV * dV;
        Udpm += kv * dV * dV;
        
        dV = kv * dV / 3.0 / V0[nc];
        for (nf = 0; nf < FNUM; nf++){
            ns1 = ns1_arr[nf];
            ns2 = ns2_arr[nf];
            ns3 = ns3_arr[nf];

            fx[ns1] += dV * R_cross_32_x[nf];
            fy[ns1] += dV * R_cross_32_y[nf];
            fz[ns1] += dV * R_cross_32_z[nf];
            fx[ns2] += dV * R_cross_13_x[nf];
            fy[ns2] += dV * R_cross_13_y[nf];
            fz[ns2] += dV * R_cross_13_z[nf];
            fx[ns3] += dV * R_cross_21_x[nf];
            fy[ns3] += dV * R_cross_21_y[nf];
            fz[ns3] += dV * R_cross_21_z[nf];
        }
        
        // Bending force
        if (kb > 0.0){
            for (nf = 0; nf < FNUM; nf++){
                A_sq = A[nf] * A[nf];
                A_rescale_x[nf] = A_vec_x[nf] / A_sq;
                A_rescale_y[nf] = A_vec_y[nf] / A_sq;
                A_rescale_z[nf] = A_vec_z[nf] / A_sq;
            }
            for (ne = 0; ne < ENUM; ne++){
                nf1 = edgelist[ne][2];
                nf2 = edgelist[ne][3];
                for (i = 0; i < 3; i++){
                    face1[i] = f_unit[nf1][i];
                    face2[i] = f_unit[nf2][i];
                }
                ns1 = edgelist[ne][0] + s_idx_start;
                ns4 = face2[3 - edgelist[ne][5]] + s_idx_start;
                x1 = px[ns1]-cx; y1 = py[ns1]-cy; z1 = pz[ns1]-cz;
                x4 = px[ns4]-cx; y4 = py[ns4]-cy; z4 = pz[ns4]-cz;
                
                y1 -= L[1] * round(y1 * L1_inv);z1 -= L[2] * round(z1 * L2_inv);
                y4 -= L[1] * round(y4 * L1_inv);z4 -= L[2] * round(z4 * L2_inv);
                
                dx14 = x4 - x1;
                dy14 = y4 - y1;
                dz14 = z4 - z1;

                //cosine
                A_dot = (A_vec_x[nf1] * A_vec_x[nf2] + A_vec_y[nf1] * A_vec_y[nf2] + A_vec_z[nf1] * A_vec_z[nf2]) / \
                        A[nf1] / A[nf2];
                if (A_dot >= 1.0){
                    theta = 0.0;
                }
                else{
                    theta = std::acos(A_dot);
                }
                if ((A_vec_x[nf1] * dx14 + A_vec_y[nf1] * dy14 + A_vec_z[nf1] * dz14) < 0.0){
                    theta = -theta;
                }

                U += kb * theta * theta;
                Udpm += kb * theta * theta;

                ns2 = face1[3 - edgelist[ne][4]] + s_idx_start;
                ns3 = edgelist[ne][1] + s_idx_start;

                x2 = px[ns2]-cx; y2 = py[ns2]-cy; z2 = pz[ns2]-cz;
                x3 = px[ns3]-cx; y3 = py[ns3]-cy; z3 = pz[ns3]-cz;

                y2 -= L[1] * round(y2 * L1_inv);z2 -= L[2] * round(z2 * L2_inv);
                y3 -= L[1] * round(y3 * L1_inv);z3 -= L[2] * round(z3 * L2_inv);
                
                dx12 = x2 - x1;
                dy12 = y2 - y1;
                dz12 = z2 - z1;
                dx13 = x3 - x1;
                dy13 = y3 - y1;
                dz13 = z3 - z1;
                dx32 = x2 - x3;
                dy32 = y2 - y3;
                dz32 = z2 - z3;
                dx34 = x4 - x3;
                dy34 = y4 - y3;
                dz34 = z4 - z3;

                r13 = std::sqrt(dx13 * dx13 + dy13 * dy13 + dz13 * dz13);
                r13_inv = kb * theta / r13;
                r13 *= (theta * kb);
                r_dot_1 = dx13 * dx32 + dy13 * dy32 + dz13 * dz32; // R13 * R32
                r_dot_2 = dx13 * dx34 + dy13 * dy34 + dz13 * dz34; // R13 * R34
                r_dot_3 = dx13 * dx12 + dy13 * dy12 + dz13 * dz12; // R13 * R12
                r_dot_4 = dx13 * dx14 + dy13 * dy14 + dz13 * dz14; // R13 * R14

                fx[ns1] -= (r_dot_1 * A_rescale_x[nf1] + r_dot_2 * A_rescale_x[nf2]) * r13_inv;
                fy[ns1] -= (r_dot_1 * A_rescale_y[nf1] + r_dot_2 * A_rescale_y[nf2]) * r13_inv;
                fz[ns1] -= (r_dot_1 * A_rescale_z[nf1] + r_dot_2 * A_rescale_z[nf2]) * r13_inv;
                fx[ns2] -= A_rescale_x[nf1] * r13;
                fy[ns2] -= A_rescale_y[nf1] * r13;
                fz[ns2] -= A_rescale_z[nf1] * r13;
                fx[ns3] += (r_dot_3 * A_rescale_x[nf1] + r_dot_4 * A_rescale_x[nf2]) * r13_inv;
                fy[ns3] += (r_dot_3 * A_rescale_y[nf1] + r_dot_4 * A_rescale_y[nf2]) * r13_inv;
                fz[ns3] += (r_dot_3 * A_rescale_z[nf1] + r_dot_4 * A_rescale_z[nf2]) * r13_inv;
                fx[ns4] -= A_rescale_x[nf2] * r13;
                fy[ns4] -= A_rescale_y[nf2] * r13;
                fz[ns4] -= A_rescale_z[nf2] * r13;
            }
        }
        //cout << E_dpm << "   ";
        /*
        s_idx_start = (nc-tN) * SNUM + tN;
        for (i = 0; i < NN_NUM; i++){
            ns = nnlist[i][0] + s_idx_start;
            ms = nnlist[i][1] + s_idx_start;
            dx13 = px[ms] - px[ns];
            if (std::abs(dx13) < D0_inter){
                dy13 = py[ms] - py[ns];
                if (std::abs(dy13) < D0_inter){
                    dz13 = pz[ms] - pz[ns];
                    if (std::abs(dz13) < D0_inter){
                        dnm = dx13 * dx13 + dy13 * dy13 + dz13 * dz13;
                        if(dnm < D0_inter * D0_inter){
                            dnm = sqrt(dnm);
                            dd = 1.0 - dnm / D0_inter;
                            //Eval += dd * dd * KC_half;  // cell-cell PE
                            ftmp = kc * (D0_inter / dnm - 1.0) / D0_inter/D0_inter;

                            dFx = ftmp * dx13;
                            dFy = ftmp * dy13;
                            dFz = ftmp * dz13;
                            fx[ns] -= dFx;
                            fy[ns] -= dFy;
                            fz[ns] -= dFz;
                            fx[ms] += dFx;
                            fy[ms] += dFy;
                            fz[ms] += dFz;
                        }
                    }
                }
            }
        }
        */
    }
    
    for (ns = tN; ns < NVTOT; ns++){
        s_idx_start = 3 * ns;
        F[s_idx_start]         += fx[ns];
        F[s_idx_start + 1]     += fy[ns];
        F[s_idx_start + 2]     += fz[ns];
    }
    
}

void tumor3D::repulsiveTumorInterfaceForces() {
    // local variables
    int ci, cj, gi, gj, bi, bj, pi, pj;
    double sij, rij, xij, dx, dy, dz, xi, ri;
    double ftmp, fx, fy, fz;
    double invRij;
    double L1_inv = 1.0 / L[1];
    double L2_inv = 1.0 / L[2];
    // loop over boxes in neighbor linked list
    for (bi = 0; bi < NBX; bi++) {
        // get start of list of vertices
        pi = head[bi];

        // loop over linked list
        while (pi > 0) {
            // real particle index
            gi = pi - 1;

            // check boundary forces
            xi = x[NDIM*gi];
            ri = r[gi];
            if (xi - wpos < ri){
                // update forces
                xij = (xi-wpos)/ri;
                fx = kc*(1.0 - xij)/ri;
                F[NDIM*gi] += fx;
                
                U += 0.5*kc*(1.0 - xij) * (1.0 - xij);
                // update wall stress
                wpress[0] += fx;
            }
            else if (xi > L[0] - ri){
                // update forces
                xij = (L[0] - xi)/ri;
                fx = -kc*(1.0 - xij)/ri;
                F[NDIM*gi] += fx;
                
                U += 0.5*kc*(1.0 - xij) * (1.0 - xij);
            }

            // cell index of gi
            ci=C_list[gi];
            
            // next particle in list
            pj = list[pi];

            // loop down neighbors of pi in same cell
            while (pj > 0) {
                // real index of pj
                gj = pj - 1;
                                
                // cell index of gj
                cj=C_list[gj];

                //do not include interaction on the same cell
                if (cj == ci) {
                    pj = list[pj];
                    continue;
                }
                
                // contact distance
                sij = r[gi] + r[gj];
                
                // particle distance
                dx = x[NDIM * gj] - x[NDIM * gi];
                if (abs(dx) < sij) {
                    dy = x[NDIM * gj + 1] - x[NDIM * gi + 1];
                    dy -= L[1] * round(dy * L1_inv);
                    if (abs(dy) < sij) {
                        dz = x[NDIM * gj + 2] - x[NDIM * gi + 2];
                        dz -= L[2] * round(dz * L2_inv);
                        if (abs(dz) < sij) {
                            rij = dx * dx + dy * dy + dz * dz;
                            if (rij < sij * sij) {
                                rij = sqrt(rij);
                                // scaled distance
                                xij = rij/sij;
                                
                                // force magnitude
                                ftmp = kc*(1 - xij)/sij;
                                
                                // increase potential energy
                                U += 0.5*kc*(1.0 - xij) * (1.0 - xij);
                                
                                // force elements
                                invRij = 1.0 / rij;
                                fx = ftmp * dx * invRij;
                                fy = ftmp * dy * invRij;
                                fz = ftmp * dz * invRij;
                                
                                // add to forces
                                F[NDIM*gi]             -= fx;
                                F[NDIM*gi + 1]         -= fy;
                                F[NDIM*gi + 2]         -= fz;
                                
                                F[NDIM*gj]             += fx;
                                F[NDIM*gj + 1]         += fy;
                                F[NDIM*gj + 2]         += fz;
                                
                            }
                        }
                    }
                }

                // update pj
                pj = list[pj];
            }

            // test overlaps with forward neighboring cells
            for (bj = 0; bj < NNN; bj++) {
                // only check if boundaries permit
                if (nn[bi][bj] == -1)
                    continue;

                // get first particle in neighboring cell
                pj = head[nn[bi][bj]];

                // loop down neighbors of pi in same cell
                while (pj > 0) {
                    // real index of pj
                    gj = pj - 1;

                    // cell index of gj
                    cj=C_list[gj];

                    //do not include interaction on the same cell
                    if (cj == ci) {
                        pj = list[pj];
                        continue;
                    }
                    
                    // contact distance
                    sij = r[gi] + r[gj];

                    // particle distance
                    dx = x[NDIM * gj] - x[NDIM * gi];
                    if (abs(dx) < sij) {
                        dy = x[NDIM * gj + 1] - x[NDIM * gi + 1];
                        dy -= L[1] * round(dy * L1_inv);
                        if (abs(dy) < sij) {
                            dz = x[NDIM * gj + 2] - x[NDIM * gi + 2];
                            dz -= L[2] * round(dz * L2_inv);
                            if (abs(dz) < sij) {
                                rij = dx * dx + dy * dy + dz * dz;
                                if (rij < sij * sij) {
                                    rij = sqrt(rij);
                                    // scaled distance
                                    xij = rij/sij;
                                    
                                    // force magnitude
                                    ftmp = kc*(1 - xij)/sij;
                                    
                                    // increase potential energy
                                    U += 0.5*kc*(1.0 - xij) * (1.0 - xij);
                                    
                                    // force elements
                                    invRij = 1.0 / rij;
                                    fx = ftmp * dx * invRij;
                                    fy = ftmp * dy * invRij;
                                    fz = ftmp * dz * invRij;
                                    
                                    // add to forces
                                    F[NDIM*gi]             -= fx;
                                    F[NDIM*gi + 1]         -= fy;
                                    F[NDIM*gi + 2]         -= fz;
                                    
                                    F[NDIM*gj]             += fx;
                                    F[NDIM*gj + 1]         += fy;
                                    F[NDIM*gj + 2]         += fz;
                                    
                                }
                            }
                        }
                    }

                    // update pj
                    pj = list[pj];
                }
            }

            // update pi index to be next
            pi = list[pi];
        }
    }
    
    wpress[0] = wpress[0] * L1_inv * L2_inv;
}

void tumor3D::stickyTumorInterfaceForces(){
    // local variables
    int ci, cj, gi, gj, bi, bj, pi, pj;
    double sij, rij, dx, dy, dz, xi, ri;
    double ftmp, fx, fy, fz;
    double invRij;
    double L1_inv = 1.0 / L[1];
    double L2_inv = 1.0 / L[2];

    // attraction shell parameters
    double shellij, cutij, xij, kint = (kc*l1)/(l2 - l1);
    
    // loop over boxes in neighbor linked list
    for (bi = 0; bi < NBX; bi++) {
        // get start of list of vertices
        pi = head[bi];
        // loop over linked list
        while (pi > 0) {
            // real particle index
            gi = pi - 1;

            // check boundary forces
            xi = x[NDIM*gi];
            ri = r[gi];
            if (xi-wpos < ri){
                // update forces
                xij = (xi-wpos)/ri;
                fx = kc*(1.0 - xij)/ri;
                F[NDIM*gi] += fx;
                
                U += 0.5*kc*(1.0 - xij) * (1.0 - xij);
                // update wall stress
                wpress[0] += fx;
            }
            else if (xi > L[0] - ri){
                // update forces
                xij = (L[0] - xi)/ri;
                fx = -kc*(1.0 - xij)/ri;
                F[NDIM*gi] += fx;
                
                U += 0.5*kc*(1.0 - xij) * (1.0 - xij);
            }
            
            // cell index of gi
            ci=C_list[gi];
            // next particle in list
            pj = list[pi];
            // loop down neighbors of pi in same cell
            while (pj > 0) {
                // real index of pj
                gj = pj - 1;

                // cell index of gj
                cj=C_list[gj];

                //do not include interaction on the same cell
                if (cj == ci) {
                    pj = list[pj];
                    continue;
                }
                
                // contact distance
                sij = r[gi] + r[gj];
                
                // attraction distances
                shellij = (1.0 + l2)*sij;
                cutij = (1.0 + l1)*sij;
                
                // particle distance
                dx = x[NDIM * gj] - x[NDIM * gi];
                if (abs(dx) < shellij) {
                    dy = x[NDIM * gj + 1] - x[NDIM * gi + 1];
                    dy -= L[1] * round(dy * L1_inv);
                    if (abs(dy) < shellij) {
                        dz = x[NDIM * gj + 2] - x[NDIM * gi + 2];
                        dz -= L[2] * round(dz * L2_inv);
                        if (abs(dz) < sij) {
                            rij = dx * dx + dy * dy + dz * dz;
                            if (rij < shellij * shellij) {
                                rij = sqrt(rij);
                                // scaled distance
                                xij = rij/sij;
                                // pick force based on vertex-vertex distance
                                // only tumor cells attract
                                if (rij > cutij && ci < tN && cj < tN){
                                    // force scale
                                    ftmp = kint*(xij - 1.0 - l2)/sij;
                                    
                                    // increase potential energy
                                    U += -0.5*kint*(1.0 + l2 - xij) * (1.0  + l2 - xij);
                                    //Utest+=-0.5*kint*(1.0 + l2 - xij) * (1.0  + l2 - xij);
                                }
                                else{
                                    if ((ci < tN && cj < tN) || rij < sij){
                                        // force scale
                                        ftmp = kc*(1 - xij)/sij;
                                        
                                        // increase potential energy
                                        if(ci < tN && cj < tN){
                                            U += 0.5*kc*((1.0 - xij) * (1.0 - xij) - l1*l2);
                                            //Utest += 0.5*kc*((1.0 - xij) * (1.0 - xij) - l1*l2);
                                        }
                                        else{
                                            U += 0.5*kc*(1.0 - xij) * (1.0 - xij);
                                        }
                                    }
                                    else{
                                        pj = list[pj];
                                        continue;
                                    }
                                }
                                // force elements
                                invRij = 1.0 / rij;
                                fx = ftmp * dx * invRij;
                                fy = ftmp * dy * invRij;
                                fz = ftmp * dz * invRij;
                                
                                // add to forces
                                F[NDIM*gi]             -= fx;
                                F[NDIM*gi + 1]         -= fy;
                                F[NDIM*gi + 2]         -= fz;
                                
                                F[NDIM*gj]             += fx;
                                F[NDIM*gj + 1]         += fy;
                                F[NDIM*gj + 2]         += fz;
                            }
                        }
                    }
                }

                // update pj
                pj = list[pj];
            }

            // test overlaps with forward neighboring cells
            for (bj = 0; bj < NNN; bj++) {
                //break;
                // only check if boundaries permit
                if (nn[bi][bj] == -1)
                    continue;

                // get first particle in neighboring cell
                pj = head[nn[bi][bj]];

                // loop down neighbors of pi in same cell
                while (pj > 0) {
                    // real index of pj
                    gj = pj - 1;

                    // cell index of gj
                    cj=C_list[gj];

                    //do not include interaction on the same cell
                    if (cj == ci) {
                        pj = list[pj];
                        continue;
                    }
                    
                    // contact distance
                    sij = r[gi] + r[gj];

                    // attraction distances
                    shellij = (1.0 + l2)*sij;
                    cutij = (1.0 + l1)*sij;

                    dx = x[NDIM * gj] - x[NDIM * gi];
                    if (abs(dx) < shellij) {
                        dy = x[NDIM * gj + 1] - x[NDIM * gi + 1];
                        dy -= L[1] * round(dy * L1_inv);
                        if (abs(dy) < shellij) {
                            dz = x[NDIM * gj + 2] - x[NDIM * gi + 2];
                            dz -= L[2] * round(dz * L2_inv);
                            if (abs(dz) < sij) {
                                rij = dx * dx + dy * dy + dz * dz;
                                if (rij < shellij * shellij) {
                                    rij = sqrt(rij);
                                    // scaled distance
                                    xij = rij/sij;
                                    
                                    // pick force based on vertex-vertex distance
                                    // only tumor cells attract
                                    if (rij > cutij && ci < tN && cj < tN){
                                        // force scale
                                        ftmp = kint*(xij - 1.0 - l2)/sij;
                                        
                                        // increase potential energy
                                        U += -0.5*kint*(1.0 + l2 - xij) * (1.0  + l2 - xij);
                                        //Utest+=-0.5*kint*(1.0 + l2 - xij) * (1.0  + l2 - xij);
                                    }
                                    else{
                                        if ((ci < tN && cj < tN) || rij < sij){
                                            // force scale
                                            ftmp = kc*(1 - xij)/sij;
                                            
                                            // increase potential energy
                                            if(ci < tN && cj < tN){
                                                U += 0.5*kc*((1.0 - xij) * (1.0 - xij) - l1*l2);
                                                //Utest += 0.5*kc*((1.0 - xij) * (1.0 - xij) - l1*l2);
                                            }
                                            else{
                                                U += 0.5*kc*(1.0 - xij) * (1.0 - xij);
                                            }
                                        }
                                        else{
                                            pj = list[pj];
                                            continue;
                                        }
                                    }
                                    // force elements
                                    invRij = 1.0 / rij;
                                    fx = ftmp * dx * invRij;
                                    fy = ftmp * dy * invRij;
                                    fz = ftmp * dz * invRij;
                                    
                                    // add to forces
                                    F[NDIM*gi]             -= fx;
                                    F[NDIM*gi + 1]         -= fy;
                                    F[NDIM*gi + 2]         -= fz;
                                    
                                    F[NDIM*gj]             += fx;
                                    F[NDIM*gj + 1]         += fy;
                                    F[NDIM*gj + 2]         += fz;
                                    
                                }
                            }
                        }
                    }

                    // update pj
                    pj = list[pj];
                }
            }

            // update pi index to be next
            pi = list[pi];
        }
    }

    wpress[0] = wpress[0] * L1_inv * L2_inv;

}

void tumor3D::repulsiveTumorInterfaceForceUpdate() {
    resetForcesAndEnergy();
    crawlerUpdate();
    repulsiveTumorInterfaceForces();
    tumorShapeForces();
}

void tumor3D::stickyTumorInterfaceForceUpdate() {
    resetForcesAndEnergy();
    crawlerUpdate();
    stickyTumorInterfaceForces();
    tumorShapeForces();
}





/******************************

    T U M O R

    F I R E

*******************************/


void tumor3D::tumorFIRE(tumor3DMemFn forceCall, double Ftol, double dt0, int div) {
    // local variables
    int i;
    // check to see if cell linked-list has been initialized
    if (NBX == -1) {
        cerr << "    ** ERROR: In dpm::fire, NBX = -1, so cell linked-list has not yet been initialized. Ending here.\n";
        exit(1);
    }

    // FIRE variables
    double P, fnorm, fcheck, vnorm, alpha, dtmax, dtmin;
    int npPos, npNeg, fireit;
    double x_max, x_min, Ldiv, xi;
    double V_wall = 0.0, F_wall = 0.0;
    //double P0=0.001;
    // set dt based on geometric parameters
    setdt(dt0);

    // Initialize FIRE variables
    alpha = alpha0;

    dtmax = 10.0 * dt;
    //1e-2
    dtmin = 1e-4 * dt;

    npPos = 0;
    npNeg = 0;

    fireit = 0;
    fcheck = 10 * Ftol;
    
    // reset forces and velocities
    resetForcesAndEnergy();
    fill(v.begin(), v.end(), 0.0);

    //Ldiv
    x_max = 0.0;
    for (i = 0; i < tN; i++){
        x_max = x_max > x[NDIM*i] ? x_max : x[NDIM*i];
    }
    x_min = L[0];
    for (i = tN; i < NVTOT; i++){
        x_min = x_min < x[NDIM*i] ? x_min : x[NDIM*i];
    }
    Ldiv = (x_max * r[tN] + x_min* r[tN-1])/(r[1] + r[tN]);
    
    // relax forces using FIRE
    while ((fcheck > Ftol || fireit < NDELAY) && fireit < itmax) {
        //printTumorInterface(0);
        // compute P
        P = 0.0;
        for (i = 0; i < vertDOF; i++)
            P += v[i] * F[i];
        P += V_wall * F_wall;
        // print to console
        
        if (fireit % NSKIP == 0) {
            cout << endl
                 << endl;
            cout << "===========================================" << endl;
            cout << "     F I R E                         " << endl;
            cout << "        M I N I M I Z A T I O N     " << endl;
            cout << "===========================================" << endl;
            cout << endl;
            cout << "    ** fireit     = " << fireit << endl;
            cout << "    ** fcheck     = " << fcheck << endl;
            cout << "   ** Ftol     = " << Ftol << endl;
            cout << "    ** Ldiv         = " << Ldiv << endl;
            cout << "    ** dt         = " << dt << endl;
            cout << "    ** P         = " << P << endl;
            cout << "    ** alpha     = " << alpha << endl;
            cout << "    ** npPos     = " << npPos << endl;
            cout << "   ** npNeg     = " << npNeg << endl;
            cout << "   ** phi      = " << vertexPreferredPackingFraction3D() << endl;
            //printTumorInterface(0);
        }
        

        // Adjust simulation based on net motion of degrees of freedom
        if (P > 0) {
            // increase positive counter
            npPos++;

            // reset negative counter
            npNeg = 0;

            // alter simulation if enough positive steps have been taken
            if (npPos > NDELAY) {
                // change time step
                if (dt * finc < dtmax)
                    dt *= finc;

                // decrease alpha
                alpha *= falpha;
            }
        }
        else {
            // reset positive counter
            npPos = 0;

            // increase negative counter
            npNeg++;

            // check if simulation is stuck
            if (npNeg > NNEGMAX) {
                cerr << "    ** ERROR: During initial FIRE minimization, P < 0 for too long, so ending." << endl;
                exit(1);
            }

            // take half step backwards, reset velocities
            for (i = 0; i < vertDOF; i++) {
                // take half step backwards
                x[i] -= (abs(dt * v[i] * 0.5)>0.01) ? 0 : dt * v[i]*0.5;

                // recenter in box
                if (x[i] > L[i % NDIM] && pbc[i % NDIM])
                    x[i] -= L[i % NDIM];
                else if (x[i] < 0.0 && pbc[i % NDIM])
                    x[i] += L[i % NDIM];
                
                // reset vertex velocities
                v[i] = 0.0;
            }
            Ldiv -=0.5 * dt * V_wall;
            V_wall = 0.0;
            
            
            // decrease time step if past initial delay
            if (fireit > NDELAY) {
                // decrease time step
                if (dt * fdec > dtmin)
                    dt *= fdec;

                // reset alpha
                alpha = alpha0;
            }
        }

        
        
        
        // VV VELOCITY UPDATE #1
        for (i = 0; i < vertDOF; i++)
            v[i] += 0.5 * dt * F[i];
        V_wall += 0.5 * dt * F_wall ;
        
        
        
        
        // compute fnorm, vnorm and P
        fnorm = 0.0;
        vnorm = 0.0;
        for (i = 0; i < vertDOF; i++) {
            fnorm += F[i] * F[i];
            vnorm += v[i] * v[i];
            
            v[i] = (v[i] > 0.001) ? 0.001 : v[i];
        }
        fnorm += F_wall * F_wall;
        vnorm += V_wall * V_wall;
        
        fnorm = sqrt(fnorm);
        vnorm = sqrt(vnorm);

        
        // update velocities (s.d. vs inertial dynamics) only if forces are acting
        
        if (fnorm > 0) {
            for (i = 0; i < vertDOF; i++)
                v[i] = (1 - alpha) * v[i] + alpha * (F[i] / fnorm) * vnorm;
            V_wall = (1 - alpha) * V_wall + alpha * (F_wall / fnorm) * vnorm;
        }

        // VV POSITION UPDATE
        for (i = 0; i < vertDOF; i++) {
            // update position
            x[i] += (abs(dt * v[i])>0.01) ? 0 : dt * v[i];

            // recenter in box
            if (x[i] > L[i % NDIM] && pbc[i % NDIM])
                x[i] -= L[i % NDIM];
            else if (x[i] < 0.0 && pbc[i % NDIM])
                x[i] += L[i % NDIM];
        }
        Ldiv += dt * V_wall;

        // update forces (function passed as argument)
        // sort particles
        reCellList3D(2.0);
        reNeighborList3D();
        CALL_MEMBER_FN(*this, forceCall)();
        F_wall=0.0;
        //Ldiv
        if(div){
            for (i = 0; i < NVTOT; i++){
                xi = x[NDIM*i];
                if (i < tN) {
                    if (xi > Ldiv - r[i]){
                        F[NDIM*i] -= kc*(1.0 - ((Ldiv - xi)/r[i]))/r[i];
                        F_wall += kc*(1.0 - ((Ldiv - xi)/r[i]))/r[i];
                    }
                        
                }
                else {
                    if (xi < Ldiv + r[i]){
                        F[NDIM*i] += kc*(1.0 - ((xi - Ldiv)/r[i]))/r[i];
                        F_wall -= kc*(1.0 - ((xi - Ldiv)/r[i]))/r[i];
                    }
                }
            }
        }
        //F_wall = (P0 - wpress[0]) * L[1] * L[2];
        // VV VELOCITY UPDATE #2
        for (i = 0; i < vertDOF; i++){
            v[i] += 0.5 * F[i] * dt;
        }
        V_wall += 0.5 * F_wall * dt;
        
        // update fcheck based on fnorm (= force per degree of freedom)
        /*
        fcheck = 0.0;
        for (i = 0; i < vertDOF; i++)
            fcheck += F[i] * F[i];
        fcheck = sqrt(fcheck / vertDOF);
         */
        fcheck = 0.0;
        for (i = 0; i < vertDOF; i++){
            if (fcheck<abs(F[i])) {
                fcheck=abs(F[i]);
            }
        }
        if (fcheck<abs(F_wall)) {
            fcheck=abs(F_wall);
        }
        
        // update iterator
        fireit++;

    }
    // check if FIRE converged
    if (fireit == itmax) {
        cout << "    ** FIRE minimization did not converge, fireit = " << fireit << ", itmax = " << itmax << "; ending." << endl;
        //warning
        //exit(1);
    }
    else {
        /*
        cout << endl;
        cout << "===========================================" << endl;
        cout << "     F I R E                         " << endl;
        cout << "        M I N I M I Z A T I O N     " << endl;
        cout << "    C O N V E R G E D!                 " << endl;
        cout << "===========================================" << endl;
        cout << endl;
        cout << "    ** fireit     = " << fireit << endl;
        cout << "    ** fcheck     = " << fcheck << endl;
        cout << "    ** U         = " << U << endl;
        cout << "    ** fnorm    = " << fnorm << endl;
        cout << "    ** vnorm     = " << vnorm << endl;
        cout << "    ** dt         = " << dt << endl;
        cout << "    ** P         = " << P << endl;
        cout << "    ** alpha     = " << alpha << endl;
        cout << endl << endl;
         */
    }
}



/******************************

    P R O T O C O L S

*******************************/


void tumor3D::setupCheck(){
    // check NVTOT set up
    if (NVTOT <= 0){
        cerr << "** ERROR: in setupCheck, NVTOT = " << NVTOT << ". Ending. " << endl;
        exit(1);
    }

    // check tN
    if (tN > NCELLS){
        cerr << "** ERROR: in setupCheck, tN = " << tN << ", which is > NCELLS = " << NCELLS << ". Ending. " << endl;
        exit(1);
    }


    if (V0.size() != NCELLS){
        cerr << "** ERROR: in setupCheck, V0.size = " << V0.size() << ", which is != NCELLS = " << NCELLS << ". Ending. " << endl;
        exit(1);
    }
    
    if (x.size() != vertDOF){
        cerr << "** ERROR: in setupCheck, V0.size = " << x.size() << ", which is != vertDOF = " << vertDOF << ". Ending. " << endl;
        exit(1);
    }
    
    if (v.size() != vertDOF){
        cerr << "** ERROR: in setupCheck, V0.size = " << v.size() << ", which is != vertDOF = " << vertDOF << ". Ending. " << endl;
        exit(1);
    }
    
    if (F.size() != vertDOF){
        cerr << "** ERROR: in setupCheck, V0.size = " << F.size() << ", which is != vertDOF = " << vertDOF << ". Ending. " << endl;
        exit(1);
    }
}


// compression with boundary forces
void tumor3D::tumorCompression(double Ftol, double Ptol, double dt0, double dphi0){
    int ci, vi, s_idx;
    double fx, fy, fz, ftmp, cx, cy, cz, dx, dy ,dz, dr_min, dr_max, r_min, r_max, R;
    vector<int> tN_list;
    int div=1;
    wpos = 0.0;
    // check correct setup
    setupCheck();

    // local variables
    int k = 0;
    double pcheck=0.0, phi0, scaleFactor = 1.0;

    //warning
    //printTumorInterface(0.0);
    //updateECMAttachments(1);
    //pcheck < Ptol
    /*
    scaleFactor=0.99;
    dr_min = 2.1 * (*max_element(r.begin(),r.end()));
    dr_max = 1.01;
    r_max=10*dr_min;
    while (r_max>dr_max) {
        for(ci = 0; ci < NCELLS; ci ++){
            if (ci<tN) {
                r_min = NN3D(ci);
            }
            else{
                r_min = 10*dr_max;
                s_idx = szList[ci];
                for(vi=0;vi<SNUM;vi++){
                    R=NN3D(s_idx + vi);
                    r_min = r_min<R ? r_min : R;
                }
            }
            while (r_min > dr_max) {
                cout <<  ci << "/" << NCELLS << ":"<<r_min << endl;
                if (ci<tN) {
                    com3D(ci, cx, cy, cz);
                    R=sqrt(cx*cx+cy*cy+cz*cz);
                    
                    if(cx < r[ci]){
                        break;
                    }
                    x[ci*NDIM] -=dr_min/R*cx*0.01;
                    x[ci*NDIM+1] -=dr_min/R*cy*0.01;
                    x[ci*NDIM+2] -=dr_min/R*cz*0.01;
                    
                    R=NN3D(ci);
                    r_min = r_min<R ? r_min : R;
                }
                else{
                    s_idx = szList[ci] * NDIM;
                    com3D(ci, cx, cy, cz);
                    if(cx < 1.0){
                        break;
                    }
                    R=sqrt(cx*cx+cy*cy+cz*cz);
                    dx = dr_min/R*cx*0.01;
                    dy = dr_min/R*cy*0.01;
                    dz = dr_min/R*cz*0.01;
                    for(vi=0;vi<SNUM;vi++){
                        x[s_idx + vi*NDIM] -=dx;
                        x[s_idx + vi*NDIM + 1] -=dy;
                        x[s_idx + vi*NDIM + 2] -=dz;
                    }
                    
                    s_idx = szList[ci];
                    for(vi=0;vi<SNUM;vi++){
                        R=NN3D(s_idx + vi);
                        r_min = r_min<R ? r_min : R;
                    }
                }
            }
        }
        r_max = 0.1*dr_max;
        for(ci = 0; ci < NCELLS; ci ++){
            if (ci<tN) {
                r_min = NN3D(ci);
            }
            else{
                r_min = 10*dr_max;
                s_idx = szList[ci];
                for(vi=0;vi<SNUM;vi++){
                    R=NN3D(s_idx + vi);
                    r_min = r_min<R ? r_min : R;
                }
            }
            r_max = r_max > r_min ? r_max : r_min;
        }
        k++;
        cout << k << ": rmax = " << r_max << endl;
        if(k>2)
            break;
    }
     */
    printTumorInterface(0.0);
    
    k=0;
    while (pcheck <Ptol && k < itmax) {
        // relax configuration (pass repsulive force update member function)
        // scale particle sizes
        // update packing fraction
        phi0 = vertexPreferredPackingFraction3D();

        if(phi0 < 0.1) {
            scaleParticleSizes3D(pow((phi0 + dphi0) / phi0,-1.0/3.0));
        }
        else if(phi0<0.35){
            L[0] *= 0.99;
            L[1] *= 0.99;
            L[2] *= 0.99;
        }
        else if(phi0<1.0){
            L[0] *= 0.999;
            L[1] *= 0.999;
            L[2] *= 0.999;
        }
        else{
            cout << " ** phi0 > 1.0. Compresssion might failed. Exit." << endl;
            exit(-1);
        }
        //else{
           // scaleFactor = pow((phi0 + dphi0) / phi0,-1.0/3.0);
            //scaleParticleSizes3D(scaleFactor);
        //}
        tumorFIRE(&tumor3D::repulsiveTumorInterfaceForceUpdate, Ftol, dt0, div);
        
        pcheck = wpress[0];
        
        // output to console
        cout << "&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&" << endl;
        cout << "===============================================" << endl << endl;
        cout << "     T U M O R    2 D                              " << endl;
        cout << "           I S O T R O P I C                         " << endl;
        cout << "            C O M P R E S S I O N                 " << endl << endl;
        cout << "===============================================" << endl;
        cout << "&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&" << endl;
        cout << endl;
        cout << "    * k             = " << k << endl;
        cout << "    * phi0             = " << phi0 << endl;
        cout << "    * scaleFactor     = " << scaleFactor << endl;
        cout << "    * pcheck         = " << pcheck << endl;
        cout << "    * U              = " << U << endl;
        cout << endl << endl;
        
    
        printTumorInterface(0.0);
        // update iterate
        k++;
    }
    
    
    
    fill(v.begin(), v.end(), 0.0);
}

// invasion at constant pressure
void tumor3D::invasionConstP(tumor3DMemFn forceCall, double M, double P0, double g0, int NT, int NPRINTSKIP){
    
    default_random_engine gen;
    normal_distribution<double> distribution(0.0, 1.0);
    
    // check correct setup
    setupCheck();
    
    // local variables
    int k, i, ci, gi, upb = 20000;
    double t = 0.0;
    double subBoxLength = 2.5;
    double x_max=0;
    double Vy = 0.0;
    double B = 0.1;
    double H=0.0;
    double M_wall = M;

    if (M!=1.0) {
        cout << "M is not 1.0. Exit" << endl;
        exit(1);
    }
    
    double F_wall = 0.0;
    double V_wall = 0.0;
    double K=0.0;
    double dt_sqr = sqrt(dt);
    double dt_2 = dt*dt;
    double dt_15 = pow(dt,1.5);
    double sqr_3 = 1.0/sqrt(3.0);
    vector<double> r1;
    vector<double> r2;
    double sg = sqrt(2*B*v0);
    vector<int> tN_list;


    for (gi=0; gi<NVTOT*NDIM+1; gi++){
        r1.push_back(0.0);
        r2.push_back(0.0);
    }

    //allocation enough memory
    x.reserve(upb);
    v.reserve(upb);
    F.reserve(upb);
    nv.reserve(upb);
    r.reserve(upb);
    V0.reserve(upb);
    Dr.reserve(upb);
    psi.reserve(upb);
    list.reserve(upb);
    szList.reserve(upb);
    pinattach.reserve(upb);
    ifbroken.reserve(upb);
    
    reCellList3D(subBoxLength);
    reNeighborList3D();
    // initial pressure
    CALL_MEMBER_FN(*this, forceCall)();
    F_wall = (P0 - wpress[0]) * L[1] * L[2];
    
    printTumorInterface(t);
    for (k=0; k<NT; k++){
        // pbcs and reset forces
        for (i=0; i<vertDOF; i++){
            // recenter in box (only if y)
            if (i % NDIM == 1){
                if (x[i] > L[1])
                    x[i] -= L[1];
                else if (x[i] < 0)
                    x[i] += L[1];
            }
        }
        
        /*******************************************************************************************************************************/
        // update positions (Velocity Verlet, OVERDAMPED) & update velocity 1st term
        for (i=0; i<vertDOF; i++){
            r1[i] = distribution(gen);
            r2[i] = distribution(gen);
            
            v[i] = v[i] + dt*0.5*F[i] -dt*0.5*B*v[i] + 0.5*dt_sqr*sg*r1[i] - 0.125*dt_2*B*(F[i]-B*v[i]) - 0.25*dt_15*B*sg*(0.5*r1[i]+sqr_3*r2[i]);
            x[i] += dt*v[i] + dt_15*sg*0.5*sqr_3*r2[i];
        }
        r1.back() = distribution(gen);
        r2.back() = distribution(gen);

        V_wall = V_wall + dt*0.5*F_wall -dt*0.5*B*V_wall + 0.5*dt_sqr*sg*r1.back() - 0.125*dt_2*B*(F_wall-B*V_wall) - 0.25*dt_15*B*sg*(0.5*r1.back()+sqr_3*r2.back());
        wpos += dt*V_wall + dt_15*sg*0.5*sqr_3*r2.back();
        
        // update psi before update force
        psiDiffusion(k);
        //psiECM();
        
        // sort particles
        reCellList3D(subBoxLength);
        reNeighborList3D();
        // update forces
        CALL_MEMBER_FN(*this, forceCall)();
        F_wall = (P0 - wpress[0]) * L[1] * L[2];
        
        // update velocity 2nd term (Velocity Verlet, OVERDAMPED)
        for (i=0; i<vertDOF; i++)
            v[i] = v[i] + 0.5*dt*F[i] - 0.5*dt*B*v[i] + 0.5*dt_sqr*sg*r1[i] - 0.125*dt_2*B*(F[i] - B*v[i]) - 0.25*dt_15*B*sg*(0.5*r1[i]+sqr_3*r2[i]);
        V_wall = V_wall + dt*0.5*F_wall -dt*0.5*B*V_wall + 0.5*dt_sqr*sg*r1.back() - 0.125*dt_2*B*(F_wall-B*V_wall) - 0.25*dt_15*B*sg*(0.5*r1.back()+sqr_3*r2.back());

        // update time
        t += dt;
        /********************************************************************************************************************************/
        
        // print message console, print position to file
        if (((k+1) % NPRINTSKIP == 0) || k==0){
            //kinetic energy
            K=0;
            for (int i = 0; i < vertDOF; i++)
                K += v[i] * v[i];
            K *= 0.5;
            
            //H
            H = U+P0*L[1]*L[2]*(L[0]-wpos) + K + V_wall*V_wall*0.5;
            
            // make sure Vy=0
            Vy=0.0;
            for (ci=0; ci<NVTOT; ci++){
                Vy += v[NDIM*ci + 1];
            }
            Vy = Vy/NVTOT;
            for (ci=0; ci<NVTOT; ci++){
                v[NDIM*ci + 1] -=Vy;
            }
            
            //find front
            x_max = 0.0;
            for (i = 0; i<tN; i++) {
                if (x[i*NDIM] > x_max) {
                    x_max = x[i*NDIM];
                }
            }
            
            cout << endl << endl;
            cout << "===========================================" << endl;
            cout << "            invading tumor cells             " << endl;
            cout << "===========================================" << endl;
            cout << endl;
            cout << "    ** k             = " << k+1 << endl;
            cout << "    ** t             = " << t << endl;
            cout << "   ** H            = " << H << endl;
            cout << "    ** P wall        = " << wpress[0] << endl;
            cout << "    ** Lx            = " << L[0]-wpos << endl;
            cout << "   ** front          = " << x_max << endl;
            cout << "   ** E              = " << U + K + V_wall*V_wall*0.5 - wpos*P0*L[1]*L[2]  << endl;
            cout << "   ** U              = " << U << endl;
            cout << "   ** Udpm           = " << Udpm << endl;
            cout << "   ** Ul             = " << Ul << endl;
            cout << "   ** Kinetic        = " << K<< endl;
            cout << "   ** Kinetic_tumor  = " << 0<< endl;
            cout << "   ** Kinetic_wall   = " << V_wall*V_wall*0.5 << endl;
            cout << "   ** potential_wall = " << -wpos*P0*L[1]*L[2] << endl;
            cout << "   ** U_tumor        = " << Utest << endl;
            cout << "\t** PRINTING POSITIONS TO FILE... " << endl;
            
            if ((k+1) % (NPRINTSKIP*10) == 0) {
                printTumorInterface(t);
                //annealing
                /*
                if ((k+1) % (NPRINTSKIP*500) == 0 && (k+1)<NT/2.0) {
                    for (int i = 0; i < vertDOF; i++)
                        v[i] *= 2.0;
                }
                else if((k+1) % (NPRINTSKIP*500) == 0 && (k+1)>NT/2.0) {
                    for (int i = 0; i < vertDOF; i++)
                        v[i] *= 0.5;
                }
                */
            }
        }
    }
}

/******************************

    P R I N T I N G

    F U N C T I O N S

*******************************/


void tumor3D::printTumorInterface(double t){
    // local variables
    int ci, vi, gi;
    double xi, yi, zi, dx, dy, dz, Lx, Ly, Lz, cx, cy ,cz;

    // check if pos object is open
    if (!posout.is_open()) {
        cerr << "** ERROR: in printConfiguration3D, posout is not open, but function call will try to use. Ending here." << endl;
        exit(1);
    }
    else
        cout << "** In printConfiguration3D, printing particle positions to file..." << endl;

    // save box sizes
    Lx = L.at(0);
    Ly = L.at(1);
    Lz = L.at(1);

    // print information starting information
    posout << setw(w) << left << "NEWFR"
           << " " << endl;
    posout << setw(w) << left << "NUMCL" << setw(w) << left << NCELLS << setw(w) << left << tN << endl;
    posout << setw(w) << left << "TSTEP" << setw(wnum) << setprecision(pnum) << left << dt << endl;
    posout << setw(w) << left << "TCURR" << setw(wnum) << setprecision(pnum) << left << t << endl;
    posout << setw(w) << left << "PACKF" << setw(wnum) << setprecision(pnum) << left << vertexPackingFraction3D() << endl;

    // print box sizes
    posout << setw(w) << left << "BOXSZ";
    posout << setw(wnum) << setprecision(pnum) << left << Lx;
    posout << setw(wnum) << setprecision(pnum) << left << Ly;
    posout << setw(wnum) << setprecision(pnum) << left << wpos;
    posout << endl;

    // print stress info
    posout << setw(w) << left << "STRSS";
    posout << setw(wnum) << setprecision(pnum) << left << 0.0;
    posout << setw(wnum) << setprecision(pnum) << left << 0.0;
    posout << setw(wnum) << setprecision(pnum) << left << 0.0;
    posout << endl;

    // print wall stress info
    posout << setw(w) << left << "WPRSS";
    posout << setw(wnum) << setprecision(pnum) << left << wpress.at(0);
    posout << setw(wnum) << setprecision(pnum) << left << wpress.at(1);
    posout << endl;

    // print coordinate for rest of the cells
    for (ci = 0; ci < NCELLS; ci++) {
        
        com3D(ci, cx, cy, cz);

        // cell information
        posout << setw(w) << left << "CINFO";
        posout << setw(w) << left << nv.at(ci);
        posout << setw(w) << left << 0;
        posout << setw(w) << left << 0;
        posout << setw(wnum) << left << V0.at(ci);
        if (ci < tN){
            posout << setw(wnum) << left << V0.at(ci);
            posout << setw(wnum) << left << 0.33;
            posout << setw(wnum) << left << 0.0;
            posout << setw(wnum) << left << Dr.at(ci);
            
        }
        else{
            posout << setw(wnum) << left << volume(ci);
            posout << setw(wnum) << left << area(ci);
            posout << setw(wnum) << left << pinpos[NDIM*(ci-tN)];
            posout << setw(wnum) << left << pinpos[NDIM*(ci-tN) + 1];
            posout << setw(w) << left << pinattach[ci-tN];
        }
        posout << endl;

        // get initial vertex positions
        gi = gindex(ci, 0);
        xi = x.at(NDIM * gi);
        yi = x.at(NDIM * gi + 1);
        zi = x.at(NDIM * gi + 2);

        // place back in box center
        yi = fmod(yi, Ly);
        zi = fmod(zi, Lz);

        posout << setw(w) << left << "VINFO";
        posout << setw(w) << left << ci;
        posout << setw(w) << left << 0;

        // output initial vertex information
        posout << setw(wnum) << setprecision(pnum) << right << xi;
        posout << setw(wnum) << setprecision(pnum) << right << yi;
        posout << setw(wnum) << setprecision(pnum) << right << zi;
        posout << setw(wnum) << setprecision(pnum) << right << r.at(gi);
        posout << setw(wnum) << setprecision(pnum) << right << v.at(NDIM * gi);
        posout << setw(wnum) << setprecision(pnum) << right << v.at(NDIM * gi + 1);
        posout << setw(wnum) << setprecision(pnum) << right << v.at(NDIM * gi + 2);
        posout << endl;

        // vertex information for next vertices
        for (vi = 1; vi < nv.at(ci); vi++) {
            // get global vertex index for next vertex
            gi++;

            xi = x.at(NDIM * gi);
            yi = x.at(NDIM * gi + 1);
            zi = x.at(NDIM * gi + 2);

            // get next vertex positions
            dx = x.at(NDIM * gi) - xi;
            xi += dx;

            dy = x.at(NDIM * gi + 1) - yi;
            dy -= Ly * round(dy / Ly);
            yi += dy;
            
            dz = x.at(NDIM * gi + 2) - zi;
            dz -= Lz * round(dz / Lz);
            zi += dz;

            // Print indexing information
            posout << setw(w) << left << "VINFO";
            posout << setw(w) << left << ci;
            posout << setw(w) << left << vi;

            // output vertex information
            posout << setw(wnum) << setprecision(pnum) << right << xi;
            posout << setw(wnum) << setprecision(pnum) << right << yi;
            posout << setw(wnum) << setprecision(pnum) << right << zi;
            posout << setw(wnum) << setprecision(pnum) << right << r.at(gi);
            posout << setw(wnum) << setprecision(pnum) << right << v.at(NDIM * gi);
            posout << setw(wnum) << setprecision(pnum) << right << v.at(NDIM * gi + 1);
            posout << setw(wnum) << setprecision(pnum) << right << v.at(NDIM * gi + 2);
            posout << endl;
        }
    }

    // print end frame
    posout << setw(w) << left << "ENDFR" << " " << endl;
}
