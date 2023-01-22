/*

    FUNCTION DEFINITIONS for tumor2D class

    Jack Treado, 06/04/21

    ** TO DO 06/24/21
    1. Add adipocyte / tumor boundary sim initialization function
    2. Make protocol that compresses ONLY so that we have premade initial conditions, don't have to spend overhead
    3. Add constant pressure via piston attached to left-side wall, box vol. changes as cells invade
        * Need to add kinetic term to stress tensor

*/


#include "tumor2D.h"
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
tumor2D::tumor2D(string &inputFileStr,int seed) : dpm(2) {
    // open file
    ifstream inputobj(inputFileStr.c_str());
    if (!inputobj.is_open()){
        cerr << "** ERROR: In tumor2D constructor, could not open file " << inputFileStr << ", ending here. " << endl;
        exit(1);
    }

    // set variables to default
    gamtt=0.0; v0=0.0; Dr0=0.0; Ds=0.0; tau=0.0; rho0=0.0; rhot=0.0; kecm=0.0; ecmbreak=0.0; pbc[0]=0; pbc[1]=1;

    // local variables
    int nvtmp, ci, vi, i;
    double val;
    double lxtmp, lytmp, wpostmp;
    double s1, s2, s3;
    double wp1, wp2;
    double a0tmp, psitmp, l0tmp, t0tmp;
    double xtmp, ytmp, rtmp, vxtmp, vytmp;
    string inputStr;

    // LINE 1: should be NEWFR
    getline(inputobj, inputStr);
    if (inputStr.compare(0,5,"NEWFR") != 0){
        cerr << "** ERROR: In tumor2D constructor, first line of input file NOT NEWFR, first line = " << inputStr << ". Ending." << endl;
        exit(1);
    }

    // read in simulation information from header
    getline(inputobj, inputStr);
    sscanf(inputStr.c_str(),"NUMCL %d %d",&NCELLS,&tN);
    //cout << "\t ** " << inputStr << endl;

    // verify input file
    if (NCELLS< 1){
        cerr << "** ERROR: in tumor2D constructor, NCELLStmp = " << NCELLS << ". Ending here." << endl;
        exit(1);
    }
    else if (tN < 1){
        cerr << "** ERROR: in tumor2D constructor, tNtmp = " << tN << ". Ending here." << endl;
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
    
    if (wpostmp!=0) {
        if (abs(log10(abs(wpostmp)))>3 || abs(wpostmp)>4) {
            wpostmp = 0;
        }
    }
    
    L.at(0) = lxtmp;
    L.at(1) = lytmp;
    wpos = wpostmp;
    // initialize stress
    getline(inputobj, inputStr);
    sscanf(inputStr.c_str(),"STRSS %lf %lf %lf",&s1,&s2,&s3);
    //cout << "\t ** " << inputStr << endl;

    stress.at(0) = s1;
    stress.at(1) = s2;
    stress.at(2) = s3;

    // initialize wall pressure
    getline(inputobj, inputStr);
    sscanf(inputStr.c_str(),"WPRSS %lf %lf",&wp1,&wp2);
    //cout << "\t ** " << inputStr << endl;

    wpress.resize(2);
    wpress.at(0) = wp1;
    wpress.at(1) = wp2;
    
    // szList and nv (keep track of global vertex indices)
    nv.resize(NCELLS);
    szList.resize(NCELLS);
    a0.resize(NCELLS);

    // initialize ecm + crawling variables
    psi.resize(tN);
    Dr.resize(tN);
    F_ij.resize(tN * NCELLS * NDIM);
    contactTime.resize(tN * (NCELLS - tN));
    
    fill(psi.begin(), psi.end(), 0.0);
    fill(Dr.begin(), Dr.end(), 0.0);
    fill(F_ij.begin(),F_ij.end(),0.0);
    fill(contactTime.begin(),contactTime.end(),0);
    
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
            sscanf(inputStr.c_str(),"CINFO %d %*d %*d %lf %*lf %*lf %lf %*f",&nvtmp,&a0tmp,&psitmp);
            psi.at(ci) = psitmp;
        }
        else
            sscanf(inputStr.c_str(),"CINFO %d %*d %*d %lf %*lf %*lf %*lf %*lf %*lf",&nvtmp,&a0tmp);

        // print to console
        //cout << "\t ** " << inputStr << endl;

        // store in vectors
        nv.at(ci) = nvtmp;
        a0.at(ci) = a0tmp;
        // update NVTOT
        NVTOT += nvtmp;

        // increment szList
        if (ci > 0)
            szList.at(ci) = szList.at(ci-1) + nv.at(ci-1);

        // loop over vertices, store coordinates
        for (vi=0; vi<nvtmp; vi++){
            // parse vertex coordinate info
            getline(inputobj, inputStr);
            sscanf(inputStr.c_str(),"VINFO %*d %*d %lf %lf %lf %lf %lf %lf %lf",&xtmp,&ytmp,&rtmp,&l0tmp,&t0tmp,&vxtmp,&vytmp);

            // push back
            x.push_back(xtmp);
            x.push_back(ytmp);
            r.push_back(rtmp);
            l0.push_back(l0tmp);
            t0.push_back(t0tmp);
            //v.push_back(vxtmp);
            //v.push_back(vytmp);
        }
        if (nv.at(ci)==1) {
            a0.at(ci) = PI*rtmp*rtmp;
        }
    }
    vertDOF = NDIM * NVTOT;
    cout << "** FINISHED LOOPING OVER DATA, CLOSING INPUT FILE OBJ\n" << endl;

    // initialize contact network (NEED TO DO HERE, NOT DONE IN DPM CONSTRUCTOR)
    cij.resize(NCELLS * (NCELLS - 1) / 2);
    for (i = 0; i < NCELLS * (NCELLS - 1) / 2; i++)
        cij.at(i) = 0;

    
    // initialize vertex indexing
    initializeVertexIndexing2D();

    // set initial l0
    setl0_init();

    //set rho0
    for(ci=tN;ci<NCELLS;ci++){
        rho0 +=a0[ci];
    }
    rho0 = rho0/(NCELLS - tN);
    //set rhot
    rhot=1;
    for(ci=0;ci<tN;ci++){
        if (rhot>a0[ci]) {
            rhot=a0[ci];
        }
    }

    // close input file object
    inputobj.close();

    // seed random number generator
    srand48(seed);
}




/*********************************

    T U M O R   C E L L

    I N I T I A L I Z A T I O N

**********************************/


// initialize single tumor cell for growth as monolayer or single crawler
void tumor2D::initializeSingleTumorCell(){
    // local variables
    int vi;
    int dtmp;

    // only proceed if tN has not been set yet, but nv etc has
    if (tN != 0){
        cout << "\t ** ERROR: in initializeSingleTumorCell, tN = " << tN << " which is not 0, cannot reset. Ending here." << endl;
        exit(1);
    }
    if (NVTOT <= 0){
        cout << "    ** ERROR: in initializeSingleTumorCell, NVTOT not assigned. Ending here." << endl;
        exit(1);
    }
    if (vertDOF <= 0){
        cout << "    ** ERROR: in initializeSingleTumorCell, vertDOF not assigned. Ending here." << endl;
        exit(1);
    }
    else if (nv.size() == 0){
        cout << "    ** ERROR: in initializeSingleTumorCell, nv vector not assigned. Ending here." << endl;
        exit(1);
    }

    // initialize tN to 1
    tN = 1;

    // create first cell centered at origin
    for (vi=0; vi<nv.at(0); vi++){
        // length from center to vertex
        dtmp = sqrt((2.0*a0.at(0))/(nv.at(0)*sin((2.0*PI)/nv.at(0))));

        // set positions
        x.at(NDIM*vi)         = dtmp*cos((2.0*PI*vi)/nv.at(0)) + 1e-2*l0[vi]*drand48();
        x.at(NDIM*vi + 1)    = dtmp*sin((2.0*PI*vi)/nv.at(0)) + 1e-2*l0[vi]*drand48();
    }
}


// set l0_init to be l0
void tumor2D::setl0_init(){
    // check that l0 set
    if (NVTOT <= 0){
        cerr << "\t ** ERROR: in setl0_init, NVTOT = " << NVTOT << ", which is <= 0. Ending. " << endl;
        exit(1);
    }
    if (l0.size() != NVTOT){
        cerr << "\t ** ERROR: in setl0_init, l0 size = " << l0.size() << ", which is != " << NVTOT << ". Ending. " << endl;
        exit(1);
    }

    // resize l0_init
    l0_init.resize(NVTOT);

    // fill with l0
    for (int gi=0; gi<NVTOT; gi++)
        l0_init.at(gi) = l0.at(gi);
}

//set t0
void tumor2D::setdt(double dt0) {
    // local variables
    int i, ci;
    double ta, tl, tb, tmin;

    rho0 = 0;
    for(ci=tN;ci<NCELLS;ci++){
        rho0 +=a0[ci];
    }
    rho0 = rho0/(NCELLS-tN);
    // set typical time scales
    ta = sqrt(a0.at(0)) / sqrt(ka);
    /*tl = (rho0 * l0.back()) / sqrt(ka * kl);
    tb = (rho0 * l0.back()) / sqrt(ka * kb);

    // set main time scale as min
    tmin = 1e8;
    if (ta < tmin)
        tmin = ta;
    if (tl < tmin)
        tmin = tl;
    if (tb < tmin)
        tmin = tb;
     
    // set dt
    dt = dt0 * tmin;
    */
    //dt = dt0 * ta;
    dt = 0.035575623680000;
}

// initialize neighbor linked list: with wall position change
void tumor2D::reNeighborLinkedList2D(double boxLengthScale) {
    // local variables
    double llscale;
    int i, d, nntmp, scx;

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
        if (d==0) {
            sb[d] = floor((L[d]-wpos) / (boxLengthScale * llscale));
            lb[d] = (L[d]-wpos) / sb[d];
        }
        else{
            sb[d] = floor(L[d] / (boxLengthScale * llscale));
            lb[d] = L[d] / sb[d];
        }
        // count total number of cells
        NBX *= sb[d];
    }
    
    // initialize list of box nearest neighbors
    scx = sb[0];
    nn.resize(NBX);

    // loop over cells, save forward neighbors for each box
    for (i = 0; i < NBX; i++) {
        // reshape entry
        nn[i].resize(NNN);

        // right neighbor (i+1)
        nn[i][0] = (i + 1) % NBX;

        // top neighbors (i,j+1), (i+1,j+1)
        if (pbc[1]){
            // (i,j+1) w/ pbc
            nn[i][1] = (i + scx) % NBX;

            // (i+1,j+1) w/ pbc
            nn[i][2] = (nn[i][1] + 1) % NBX;
        }
        else {
            // if on top row, both = -1
            if (i >= NBX - scx){
                nn[i][1] = -1;
                nn[i][2] = -1;
            }
            // if not on top row, still add
            else{
                nn[i][1] = i + scx;
                nn[i][2] = nn[i][1] + 1;
            }
        }

        // bottom neighbor w/ pbc (j-1)
        nntmp = (i + NBX - scx) % NBX;

        // bottom-right neighbor (i+1, j-1)
        if (pbc[1])
            nn[i][3] = nntmp + 1;
        else{
            // if on bottom row, skip
            if (i < scx)
                nn[i][3] = -1;
            // otherwise, set
            else
                nn[i][3] = nntmp + 1;
        }

        // right-hand bc (periodic)
        if ((i + 1) % scx == 0) {
            if (pbc[0]) {
                nn[i][0] = i - scx + 1;
                if (pbc[1]){
                    nn[i][2] = nn[i][1] - scx + 1;
                    nn[i][3] = nntmp - scx + 1;
                }
            }
            else {
                nn[i][0] = -1;
                nn[i][2] = -1;
                nn[i][3] = -1;
            }
        }
    }

    // linked-list variables
    head.resize(NBX);
    last.resize(NBX);
    list.resize(NVTOT + 1);

    // print box info to console
    //cout << "initially NBX = " << NBX << " ..." << endl;
    //cout << "sb = " << sb[0] << "," << sb[1] << "...lb = " << lb[0] << "," << lb[1] << endl;
}

//get head, last, and list with wall position changing
void tumor2D::neighborLinkedList2D() {
    // local variables
    int d, gi, boxid, sbtmp;
    double xtmp;

    // reset linked list info
    fill(list.begin(), list.end(), 0);
    fill(head.begin(), head.end(), 0);
    fill(last.begin(), last.end(), 0);

    // sort vertices into linked list
    for (gi = 0; gi < NVTOT; gi++) {
        // 1. get cell id of current particle position
        boxid = 0;
        sbtmp = 1;
        for (d = 0; d < NDIM; d++) {
            // current location
            xtmp = x[NDIM * gi + d];

            // check out-of-bounds
            if (xtmp < 0){
                if (pbc[d])
                    xtmp -= L[d]*floor(xtmp/L[d]);
                else
                    if (d==0) {
                        if (xtmp < wpos) {
                            xtmp = wpos + 0.00001;
                        }
                    } else {
                        xtmp = 0.00001;
                    }
            }
            else if (xtmp > L[d]){
                if (pbc[d])
                    xtmp -= L[d]*floor(xtmp/L[d]);
                else
                    xtmp = 0.99999*L[d];
            }

            // add d index to 1d list
            if (d==0) {
                boxid += floor((xtmp-wpos) / lb[d]) * sbtmp;
            } else {
                boxid += floor(xtmp / lb[d]) * sbtmp;
            }
            
            // increment dimensional factor
            sbtmp *= sb[d];
        }
        //boxid = 0;
        // 2. add to head list or link within list
        // NOTE: particle ids are labelled starting from 1, setting to 0 means end of linked list
        if (head[boxid] == 0) {
            head[boxid] = gi + 1;
            last[boxid] = gi + 1;
        }
        else {
            list[last[boxid]] = gi + 1;
            last[boxid] = gi + 1;
        }
    }
}

// set monolayer positions to center of box
void tumor2D::initializeTumorMonolayerPositions(double phi0, double Ftol, double kwell){
    // local variables
    int i, d, ci, cj, vi, vj, gi, cellDOF = NDIM * NCELLS;
    double areaSum, dnorm, xtra = 1.1;

    // initialize ecm + crawling variables
    psi.resize(tN);
    Dr.resize(tN);
    F_ij.resize(tN * NCELLS * NDIM);
    contactTime.resize(tN * (NCELLS - tN));
    
    setl0_init();

    fill(psi.begin(), psi.end(), 2.0*PI*drand48());
    fill(Dr.begin(), Dr.end(), Dr0);
    fill(F_ij.begin(),F_ij.end(),0.0);
    fill(contactTime.begin(),contactTime.end(),0);
    
    // local disk vectors
    vector<double> drad(NCELLS, 0.0);
    vector<double> dpos(cellDOF, 0.0);
    vector<double> dv(cellDOF, 0.0);
    vector<double> dF(cellDOF, 0.0);

    // print to console
    cout << "** initializing particle positions using 2D SP model and FIRE relaxation ..." << endl;

    // initialize box size based on packing fraction
    areaSum = 0.0;
    for (ci = 0; ci < NCELLS; ci++){
        if (nv.at(ci)==1) {
            areaSum += a0[ci];
        } else {
            areaSum += a0[ci] + 0.25 * PI * pow(2.0 * r.at(szList[ci]), 2.0) * (0.5 * nv.at(ci) - 1);
        }
    }
    
    // set box size
    for (d = 0; d < NDIM; d++)
        L.at(d) = pow(areaSum / phi0, 1.0 / NDIM);

    // initialize cell centers randomly
    for (ci = 0; ci < cellDOF; ci += 2)
        dpos.at(ci) = L[ci % 2] * drand48();
    for (ci = cellDOF - 1; ci > 0; ci -= 2)
        dpos.at(ci) = L[ci % 2] * drand48();

    // set radii of SP disks
    for (ci = 0; ci < NCELLS; ci++)
        drad.at(ci) = xtra * sqrt((2.0 * a0.at(ci)) / (nv.at(ci) * sin(2.0 * PI / nv.at(ci))));

    // FIRE VARIABLES
    double P = 0;
    double fnorm = 0;
    double vnorm = 0;
    double alpha = alpha0;

    double dt0 = 1e-2;
    double dtmax = 10 * dt0;
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

            // also add well force toward center of box
            dnorm = sqrt(pow(dpos[NDIM*ci] - 0.5*L[0],2.0) + pow(dpos[NDIM*ci + 1] - 0.5*L[1],2.0));
            dF[NDIM * ci] -= kwell*(dpos[NDIM*ci] - 0.5*L[0])/dnorm;
            dF[NDIM * ci + 1] -= kwell*(dpos[NDIM*ci + 1] - 0.5*L[1])/dnorm;
        }

        // FIRE step 5. Final VV update
        for (i = 0; i < cellDOF; i++)
            dv[i] += 0.5 * dt * dF[i];

        // update forces to check
        fcheck = 0.0;
        for (i = 0; i < cellDOF; i++)
            fcheck += dF[i] * dF[i];
        fcheck = sqrt(fcheck / NCELLS);

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
        cout << "    C O N V E R G E D!                 " << endl<< endl;

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
    for (ci = 0; ci < NCELLS; ci++) {
        for (vi = 0; vi < nv.at(ci); vi++) {
            // get global vertex index
            gi = gindex(ci, vi);

            // length from center to vertex
            dtmp = sqrt((2.0 * a0.at(ci)) / (nv.at(ci) * sin((2.0 * PI) / nv.at(ci))));

            // set positions
            x.at(NDIM * gi) = dtmp * cos((2.0 * PI * vi) / nv.at(ci)) + dpos.at(NDIM * ci) + 1e-2 * l0[gi] * drand48();
            x.at(NDIM * gi + 1) = dtmp * sin((2.0 * PI * vi) / nv.at(ci)) + dpos.at(NDIM * ci + 1) + 1e-2 * l0[gi] * drand48();
        }
    }
}


// Initialize collection of tumor cells and adipocytes
//
// NOTE: parameters for adipocytes(tumors)
// start with a(t) in camelCase variables
//
// parameters
// -- calA0:         shape parameter
// -- disp:         size dispersity
// -- areaRatio:     ratio of adipocyte preferred areas to tumor preferred areas
// -- NV:             number of vertices on particle with average sqrt(a0)
//
void tumor2D::initializeTumorInterface(double aCalA0, double tCalA0, double aDisp, double tDisp, double areaRatio, int aNV, int tNV){
    // local variables
    int ci, nvtmp;
    double lenscale, r1, r2, grv;

    // print to console
    cout << "** initializing gaussian tumor and adipocyte DPM particles in 2D with size dispersions:" << endl;
    cout << "\t** aDisp = " << aDisp << endl;
    cout << "\t** tDisp = " << tDisp << endl;
    cout << "** setting up nv + szList, setting shape parameters and initializing indexing ..." << endl;

    // szList and nv (keep track of global vertex indices)
    nv.resize(NCELLS);
    szList.resize(NCELLS);

    // initialize ecm + crawling variables
    psi.resize(tN);
    Dr.resize(tN);
    F_ij.resize(tN * NCELLS * NDIM);
    contactTime.resize(tN * (NCELLS - tN));
    
    fill(psi.begin(), psi.end(), 0.0);
    fill(Dr.begin(), Dr.end(), 0.0);
    fill(F_ij.begin(),F_ij.end(),0.0);
    fill(contactTime.begin(),contactTime.end(),0);
    
    pinpos.resize(NDIM * (NCELLS - tN));
    pinattach.resize(NCELLS - tN);
    ifbroken.resize(NCELLS - tN);

    // initialize number of vertices on each cell
    nv.at(0) = tNV;
    NVTOT = tNV;
    for (ci=1; ci<NCELLS; ci++){
        // use Box-Muller to generate polydisperse sample
        r1 = drand48();
        r2 = drand48();
        grv = sqrt(-2.0*log(r1))*cos(2.0*PI*r2);
        if (ci < tN)
            nvtmp = tNV;//floor(tDisp*tNV*grv + tNV);
        else
            nvtmp = floor(aDisp*aNV*grv + aNV);

        //if (nvtmp < nvmin)
            //nvtmp = nvmin;

        // store size of cell ci
        nv.at(ci) = nvtmp;
        szList.at(ci) = szList.at(ci-1) + nv.at(ci-1);

        // add to total NV count
        NVTOT += nvtmp;
    }
    vertDOF = NDIM * NVTOT;

    // resize shape paramters
    l0.resize(NVTOT);
    l0_init.resize(NVTOT);
    t0.resize(NVTOT);
    r.resize(NVTOT);

    // initialize particle sizes based on areaRatio
    for (ci=0; ci<NCELLS; ci++){
        if (ci < tN){
            lenscale = (double) nv.at(ci) / tNV;
            
            
            //bidisperse
            if (ci < tN/2) {
                lenscale = lenscale * 0.9;
            }
            else{
                lenscale = lenscale * 1.1;
            }
            
                        
            initializeVertexShapeParameters(ci,tCalA0,lenscale);
        }
        else{
            lenscale = (double) (nv.at(ci) * sqrt(areaRatio)) / aNV;
            initializeVertexShapeParameters(ci,aCalA0,lenscale);
        }
    }

    // initialize l0_init
    setl0_init();

    // initialize vertex indexing
    initializeVertexIndexing2D();
}


void tumor2D::initializeTumorInterfacePositions(double phi0, double Ftol, double prt, double aspectRatio){
    // local variables
    int i, d, ci, cj, vi, vj, gi, cellDOF = NDIM * NCELLS;
    double areaSum, xtra = 1.05, xi, Ldiv;
    double aspect_ratio = aspectRatio;
    // local disk vectors
    vector<double> drad(NCELLS, 0.0);
    vector<double> dpos(cellDOF, 0.0);
    vector<double> dv(cellDOF, 0.0);
    vector<double> dF(cellDOF, 0.0);

    // print to console
    cout << "** initializing particle positions using 2D SP model and FIRE relaxation ..." << endl;

    // initialize box size based on packing fraction
    //if N=1
    areaSum = 0.0;
    for (ci = 0; ci < NCELLS; ci++){
        if (nv.at(ci)==1) {
            areaSum += a0[ci];
        } else {
            areaSum += a0[ci] + 0.25 * PI * pow(2.0 * r.at(szList[ci]), 2.0) * (0.5 * nv.at(ci) - 1);
        }
    }

    // set box size
    L.at(1) = sqrt(areaSum/(aspect_ratio*phi0));
    L.at(0) = aspect_ratio*L[1];

    // dividing wall position between adipocytes and
    Ldiv = prt * L[0];

    // initialize tumor cell centers in left-hand partition of the box
    for (ci=0; ci<tN; ci++){
        dpos.at(NDIM*ci)         = (Ldiv - 2.0*drad[ci])*drand48() + drad[ci];
        dpos.at(NDIM*ci + 1)     = L[1]*drand48();
    }

    // initialize WAT cell centers to the right
    for (ci=tN; ci<NCELLS; ci++){
        dpos.at(NDIM*ci)         = (L[0] - Ldiv - 2.0*drad[ci])*drand48() + Ldiv + drad[ci];
        dpos.at(NDIM*ci + 1)     = L[1]*drand48();
    }

    // set radii of SP disks
    //if N=1
    for (ci = 0; ci < NCELLS; ci++){
        if (nv.at(ci)==1) {
            drad.at(ci) += sqrt(a0.at(ci)/PI);
        } else {
            drad.at(ci) = xtra * sqrt((2.0 * a0.at(ci)) / (nv.at(ci) * sin(2.0 * PI / nv.at(ci))));
        }
    }

    // FIRE VARIABLES
    double P = 0;
    double fnorm = 0;
    double vnorm = 0;
    double alpha = alpha0;

    double dt0 = 1e-2;
    double dtmax = 10 * dt0;
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
                else if (xi > Ldiv - drad[ci])
                    dF[NDIM*ci] -= (1.0 - ((Ldiv - xi)/drad[ci]))/drad[ci];
            }
            else {
                if (xi < Ldiv + drad[ci])
                    dF[NDIM*ci] += (1.0 - ((xi - Ldiv)/drad[ci]))/drad[ci];
                else if (xi > L[0] - drad[ci])
                    dF[NDIM*ci] -= (1.0 - ((L[0] - xi)/drad[ci]))/drad[ci];
            }
        }

        // FIRE step 5. Final VV update
        for (i = 0; i < cellDOF; i++)
            dv[i] += 0.5 * dt * dF[i];

        // update forces to check
        fcheck = 0.0;
        for (i = 0; i < cellDOF; i++)
            fcheck += dF[i] * dF[i];
        fcheck = sqrt(fcheck / NCELLS);

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
            } else {
                
                dtmp = sqrt((2.0 * a0.at(ci)) / (nv.at(ci) * sin((2.0 * PI) / nv.at(ci))));
                
                // set positions
                x.at(NDIM * gi) = dtmp * cos((2.0 * PI * vi) / nv.at(ci)) + dpos.at(NDIM * ci) + 1e-2 * l0[gi] * drand48();
                x.at(NDIM * gi + 1) = dtmp * sin((2.0 * PI * vi) / nv.at(ci)) + dpos.at(NDIM * ci + 1) + 1e-2 * l0[gi] * drand48();
            }
        }
    }
}


/*********************************

    E D I T I N G   &

            U P D A T I N G

**********************************/

// divide a single cell, assume preallocation
void tumor2D::divide(int ci){
    // local variable
    int gi, g1tmp, g2tmp, vi, ci1, ci2, nv1, nv2, icut1, icut2, nh1, nh2, vitmp;
    double Dx, Dy, Dnorm, dhatx, dhaty, nhatx, nhaty, xtmp, ytmp;

    // cell indices
    ci1 = tN - 1;
    ci2 = ci1 + 1;

    // number of vertices on each cell
    nv1     = nv.at(ci1);
    nv2     = nv.at(ci2);

    // end point of stitch
    icut1     = szList.at(ci1) + floor(0.5*nv1) - 1;
    icut2     = szList.at(ci1) + nv1 - 1;

    // get size of semicircles
    nh1     = floor(0.5*nv1);
    nh2     = nv1 - nh1;

    // move second-half of vertices from 1 -> first half of 2
    for (gi=0; gi<nh2; gi++){
        // get temporary global vertex indices
        g1tmp = icut1 + gi + 1;
        g2tmp = szList.at(ci2) + gi;

        // make old vertices part of cell 2
        x.at(NDIM*g2tmp) = x.at(NDIM*g1tmp);
        x.at(NDIM*g2tmp + 1) = x.at(NDIM*g1tmp + 1);
    }

    // get centerline
    Dx = x.at(NDIM*szList.at(ci1)) - x.at(NDIM*icut1);
    if (pbc[0])
        Dx -= L[0]*round(Dx/L[0]);
    Dy = x.at(NDIM*szList.at(ci1) + 1) - x.at(NDIM*icut1 + 1);
    if (pbc[1])
        Dy -= L[1]*round(Dy/L[1]);
    Dnorm = sqrt(Dx*Dx + Dy*Dy);

    // centerline unit vector
    dhatx = Dx/Dnorm;
    dhaty = Dy/Dnorm;

    // centerline normal vector
    nhatx = -dhaty;
    nhaty = dhatx;

    // stitch vertices on centerline (cell 1)
    vitmp = 1;
    // for (gi=nh1+1; gi<nv1; gi++){
    //     // get new position
    //     xtmp +=
    // }
}







/******************************

    B I O L O G I C A L

    F U N C T I O N S

*******************************/


// -- CRAWLING CELLS

// update psi based on Dr only
void tumor2D::psiDiffusion(){
    // local variables
    int ci;
    double r1, r2, grv;

    // update director for each cell
    for (ci=0; ci<tN; ci++){
        // generate random variable
        r1 = drand48();
        r2 = drand48();
        grv = sqrt(-2.0*log(r1))*cos(2.0*PI*r2);

        // update director for cell ci
        psi[ci] += sqrt(2.0*dt*Dr0)*grv;
    }
}

//update psi based on ECM
void tumor2D::psiECM(){
    //local variables
    int ci, cj, ck, cj_max;
    double ftmp, fxtmp, fytmp, Ftangent, dpsi, sign_teller, f_max;
    double rnd_d;
    //loop over all tumor cells
    for (ci=0; ci<tN; ci++){
        
        //find max pressure
        cj_max = tN;
        f_max = 0;
        for (cj=tN; cj<NCELLS; cj++){
            if (F_ij[NDIM*(ci*NCELLS+cj)]!=0) {
                ftmp = sqrt(F_ij[NDIM*(ci*NCELLS+cj)] * F_ij[NDIM*(ci*NCELLS+cj)] + F_ij[NDIM*(ci*NCELLS+cj)+1] * F_ij[NDIM*(ci*NCELLS+cj)+1]);
                if (ftmp>f_max) {
                    f_max = ftmp;
                    cj_max = cj;
                }
            }
        }
        
        if (f_max > 0) {
            //tangential direction
            fxtmp = F_ij[NDIM*(ci*NCELLS+cj_max) + 1]/f_max;
            fytmp = -F_ij[NDIM*(ci*NCELLS+cj_max)]/f_max;
            
            //choose the side along with previous psi
            sign_teller = fxtmp * cos(psi[ci]) + fytmp * sin(psi[ci]);
            if (sign_teller < 0) {
                fxtmp = -fxtmp;
                fytmp = -fytmp;
            }
            
            //choose the side with less pressure from other tumor cells, if there any
            sign_teller = 0;
            for (ck=0; ck<tN; ck++){
                if (ck!=ci) {
                    if (F_ij[NDIM*(ci*NCELLS+ck)]!=0) {
                        sign_teller += F_ij[NDIM*(ci*NCELLS+ck)] * fxtmp + F_ij[NDIM*(ci*NCELLS+ck) + 1] * fytmp;
                    }
                }
            }
            if (sign_teller < 0) {
                fxtmp = -fxtmp;
                fytmp = -fytmp;
            }

            
            //compute angle
            /*
            Ftangent = atan2(fytmp, fxtmp);
            dpsi = Ftangent - psi[ci];
            dpsi -= 2.0*PI*round(dpsi/(2.0*PI));
            psi[ci] += dpsi;
            */
            psi[ci] = atan2(fytmp, fxtmp);
            
        }
        
    }
}

//tumor cells swim
void tumor2D::crawlerUpdate(){
    // local variables
    int gi, ci, cj, vi;
    double cx, cy, rix, riy, ux, uy, psitmp, dpsi, v0tmp, rnorm;
    // loop over all cells, vertices
    gi = 0;
    for (ci=0; ci<tN; ci++){
        // loop over vertices
        for (vi=0; vi<nv[ci]; vi++){
            F[NDIM*gi] += v0*cos(psi[ci]);
            F[NDIM*gi + 1] += v0*sin(psi[ci]);
            gi++;
        }
    }
}

//Adipocyte shrink under pressure
void tumor2D::adipocyteShrink(){
    int ci, gi, vi;
    double atmp, a0tmp, da;
    for (ci=tN; ci<NCELLS; ci++){
        atmp = area(ci);
        a0tmp = a0[ci];
        da = (atmp / a0tmp) - 1.0;
        if (da < 0 && ka * da * da > 0.01) {
            a0[ci] = atmp;
            gi = szList.at(ci);
            for (vi = 0; vi < nv[ci]; vi++) {
                l0.at(gi + vi) = l0.at(gi + vi) * sqrt(1+da);
                r.at(gi + vi) = 0.5 * l0.at(gi + vi);
            }
        }
    }
}

//tumor cell growth when in contact with adipocytes
//then adipocyte shrink too maintain costant volume
void tumor2D::tumorGrowth(double g0){
    int ci, cj, gi, gj, vi, vj;
    double a0tmp, sq_a_ratio;
    //warning
    double growth_rate = g0;
    //loop over tumor cells
    //cout << *max_element(contactTime.begin(),contactTime.end()) << endl;
    for (ci=0; ci<tN; ci++){
        //Is this tumor cell in contact with any adipocyte? loop over adipocytes
        for (cj=tN; cj<NCELLS; cj++){
            //warning
            //ftmp = 1; 0.015
            if(contactTime[ci * (NCELLS - tN) + cj - tN]>1000){
                contactTime[ci * (NCELLS - tN) + cj - tN] = 0;
                //grow
                a0[ci] += growth_rate;
                gi = szList.at(ci);
                for (vi = 0; vi < nv[ci]; vi++) {
                    r.at(gi + vi) = sqrt(a0.at(gi + vi)/PI);
                }
                //shrink
                a0[cj] -= growth_rate;
                sq_a_ratio = sqrt(a0[cj]/(a0[cj]+growth_rate));
                gj = szList.at(cj);
                for (vj = 0; vj < nv[cj]; vj++) {
                    l0.at(gj + vj) = l0.at(gj + vj) * sq_a_ratio;
                    r.at(gj + vj) = 0.5 * l0.at(gj + vj);
                }
            }
        }
    }
}

//tumor cell divide when too large
void tumor2D::tumorDivide(double g0){
    //local variables
    int ci, cj, gi, vi, vim1, vip1;
    double rnd_t, x2, y2;
    double radius_0 = sqrt(rhot/PI);
    double delta_r = (1+sqrt(2))/3*radius_0;
    
    //loop over tumor cells
    for (ci=0; ci<tN; ci++){
        //if too large, then divide
        //a0[ci]> 2*rhot
        if (drand48()<g0) {
            //put a new disk on top of it, add noise on COM of both disk, update radii
            //random
            rnd_t = 2.0*PI*drand48();

            x[NDIM*ci]     += delta_r * cos(rnd_t);
            x[NDIM*ci + 1] += delta_r * sin(rnd_t);
            
            r[ci] = radius_0;
            a0[ci] = rhot;

            x2 = x[NDIM*ci] - 2 * delta_r * cos(rnd_t);
            y2 = x[NDIM*ci+1] -2* delta_r * sin(rnd_t);
            
            //insert to the front of the original disk
            x.insert(x.begin() + NDIM*ci    , x2);
            x.insert(x.begin() + NDIM*ci + 1, y2);
            
            v.insert(v.begin() + NDIM*ci    , 0);
            v.insert(v.begin() + NDIM*ci + 1, 0);

            F.insert(F.begin() + NDIM*ci    , 0);
            F.insert(F.begin() + NDIM*ci + 1, 0);
            
            nv.insert(nv.begin() + ci, 1);
            r.insert(r.begin() + ci, radius_0*1.1/0.9);
            a0.insert(a0.begin() + ci, rhot);
            l0.insert(l0.begin() + ci, l0[0]);
            t0.insert(t0.begin() + ci, t0[0]);
            l0_init.insert(l0_init.begin() + ci, l0[0]);

            
            Dr.insert(Dr.begin() + ci, Dr0);
            psi.insert(psi.begin() + ci, 2.0*PI*drand48());
            
            //contactTime
            for (cj = 0; cj< (NCELLS - tN); cj++) {
                contactTime.insert(contactTime.begin() + (ci + 1) * (NCELLS - tN), 0);
                contactTime[ci * (NCELLS - tN) + cj] = 0;
            }
            
            //update global variables
            NCELLS += 1;
            NVTOT  += 1;
            tN     += 1;
            vertDOF = NDIM * NVTOT;
            

            //update global vectors
            fill(cij.begin(),cij.end(),0);
            fill(F_ij.begin(),F_ij.end(),0.0);
            fill(im1.begin(),im1.end(),0);
            fill(ip1.begin(),ip1.end(),0);
            fill(head.begin(),head.end(),0);
            fill(last.begin(),last.end(),0);
            fill(list.begin(),list.end(),0);
            
            cij.resize(NCELLS * (NCELLS - 1) / 2);
            F_ij.resize(tN * NCELLS * NDIM);
            head.resize(NBX);
            last.resize(NBX);
            
            
            list.push_back(0);
            //im1.resize(NVTOT);
            //ip1.resize(NVTOT);

            //list: box list, cell list
            //box list
            sortNeighborLinkedList2D();

            szList.resize(NCELLS);
            fill(szList.begin(), szList.end(), 0.0);
            for (cj = 1; cj < NCELLS; cj++) {
                szList.at(cj) = szList.at(cj-1) + nv.at(cj-1);
            }
            
            // initializeVertexIndexing2D
            im1.resize(NVTOT);
            ip1.resize(NVTOT);
            for (cj = 0; cj < NCELLS; cj++) {
                // vertex indexing
                for (vi = 0; vi < nv.at(cj); vi++) {
                    // wrap local indices
                    vim1 = (vi - 1 + nv.at(cj)) % nv.at(cj);
                    vip1 = (vi + 1) % nv.at(cj);

                    // get global wrapped indices
                    gi = gindex(cj, vi);
                    im1.at(gi) = gindex(cj, vim1);
                    ip1.at(gi) = gindex(cj, vip1);
                }
            }

            
            //warning
            break;
            
        }
    }
    //loop over adiposites
    /*
    for (ci=tN; ci<NCELLS; ci++){
        //if adipocytes are too small, then become tumor cell
        //warning
        if (a0[ci]< 2*rhot) {
            //store position before delete
            com2D(ci,x2,y2);
            
            //delete adipocyte
            gi = szList.at(ci);
            
            x.erase(x.begin() + NDIM*gi, x.begin() + NDIM*(gi+nv[ci]));
            v.erase(v.begin() + NDIM*gi, v.begin() + NDIM*(gi+nv[ci]));
            F.erase(F.begin() + NDIM*gi, F.begin() + NDIM*(gi+nv[ci]));
            
            r.erase(r.begin() + gi, r.begin() + gi+nv[ci]);
            l0.erase(l0.begin() + gi, l0.begin() + gi+nv[ci]);
            t0.erase(t0.begin() + gi, t0.begin() + gi+nv[ci]);
            l0_init.erase(l0_init.begin() + gi, l0_init.begin() + gi+nv[ci]);
            a0.erase(a0.begin() + ci);
            
            //update NV
            NVTOT  -= nv[ci];
            NVTOT  += 1;
            nv.erase(nv.begin() + ci);
            
            pinpos.erase(pinpos.begin() + NDIM*(ci - tN));
            pinpos.erase(pinpos.begin() + NDIM*(ci - tN));
            //ifbroken.erase(ifbroken.begin() + ci - tN);
            //pinattach.erase(pinattach.begin() + ci - tN);
            
            //insert tumor disk to the end of the tumor list
            x.insert(x.begin() + NDIM*(tN)    , x2);
            x.insert(x.begin() + NDIM*(tN) + 1, y2);
            
            v.insert(v.begin() + NDIM*(tN)   , 0);
            v.insert(v.begin() + NDIM*(tN) + 1, 0);

            F.insert(F.begin() + NDIM*(tN)    , 0);
            F.insert(F.begin() + NDIM*(tN) + 1, 0);
            
            nv.insert(nv.begin() + (tN), 1);
            r.insert(r.begin() + (tN), sqrt(2)*radius_0);
            a0.insert(a0.begin() + (tN), 2*rhot);
            l0.insert(l0.begin() + (tN), l0[0]);
            t0.insert(t0.begin() + (tN), t0[0]);
            l0_init.insert(l0_init.begin() + (tN), l0[0]);

            
            Dr.insert(Dr.begin() + (tN), Dr0);
            psi.insert(psi.begin() + (tN), 2.0*PI*drand48());
            
            //contactTime
            for (cj = 0; cj< (NCELLS - tN); cj++) {
                contactTime.insert(contactTime.end(), 0);
            }
            
            //use tN insdead of tN-1 since the previous loop inserted a tumor
            //erase backwards
            for (cj = tN; cj> -1; cj--) {
                contactTime.erase(contactTime.begin() + cj * (NCELLS - tN) + ci - tN);
            }

            //update global variables
            tN     += 1;
            vertDOF = NDIM * NVTOT;

            //update global vectors
            fill(cij.begin(),cij.end(),0);
            fill(F_ij.begin(),F_ij.end(),0.0);
            fill(im1.begin(),im1.end(),0);
            fill(ip1.begin(),ip1.end(),0);
            fill(head.begin(),head.end(),0);
            fill(last.begin(),last.end(),0);
            fill(list.begin(),list.end(),0);
            
            cij.resize(NCELLS * (NCELLS - 1) / 2);
            F_ij.resize(tN * NCELLS * NDIM);
            head.resize(NBX);
            last.resize(NBX);
            
            list.push_back(0);
            //im1.resize(NVTOT);
            //ip1.resize(NVTOT);

            //list: box list, cell list
            //box list
            sortNeighborLinkedList2D();

            szList.resize(NCELLS);
            fill(szList.begin(), szList.end(), 0.0);
            for (cj = 1; cj < NCELLS; cj++) {
                szList.at(cj) = szList.at(cj-1) + nv.at(cj-1);
            }
            
            //initializeVertexIndexing2D(): see tumor divide
            break;
        }
    }
    */
}

// -- ADIPOCYTE ECM ATTACHEMENT

// update
void tumor2D::updateECMAttachments(bool attach){
    // check that vectors have been created
    if (pinattach.size() != NCELLS-tN){
        cerr << "    ** ERROR: in updatePinAttachments, number of pins != number of adipocytes, so ending here." << endl;
        exit(1);
    }
    if (pinpos.size() != NDIM*(NCELLS - tN)){
        cerr << "    ** ERROR: in updatePinAttachments, number of pins != number of adipocytes, so ending here." << endl;
        exit(1);
    }

    // local variables
    int ci;
    double cx, cy;

    // set all attached
    fill(pinattach.begin(), pinattach.end(), attach);

    // update positions
    for (ci=tN; ci<NCELLS; ci++){
        com2D(ci,cx,cy);
        pinpos.at(NDIM*(ci-tN)) = cx;
        pinpos.at(NDIM*(ci-tN) + 1) = cy;
    }
}

// add to force
void tumor2D::adipocyteECMAdhesionForces(){
    int gi, ci, vi, xind, yind, nvtmp, pinCellInd;
    double dpinx, dpiny, dpin, cx, cy;

    // loop over cells
    for (ci=tN; ci<NCELLS; ci++){
        // get pin indices
        pinCellInd = ci - tN;

        // tmp number of vertices
        nvtmp = nv[ci];

        // get distance to pin if attached
        if (pinattach[pinCellInd]){
            // get centers of mass
            com2D(ci,cx,cy);

            // get distance to pin
            dpinx = pinpos[NDIM*pinCellInd] - cx;
            if (pbc[0])
                dpinx -= L[0]*round(dpinx/L[0]);

            dpiny = pinpos[NDIM*pinCellInd + 1] - cy;
            if (pbc[1])
                dpiny -= L[1]*round(dpiny/L[1]);

            dpin = sqrt(pow(dpinx,2.0) + pow(dpiny,2.0));

            // pin forces on each vertex
            gi = gindex(ci,0);
            for (vi=0; vi<nvtmp; vi++){
                // positions using global indexing
                xind = NDIM*(gi+vi);
                yind = xind + 1;

                // if an adipocyte and pin is intact, compute force due to pinning spring
                if (dpin < ecmbreak*sqrt(a0[ci])){
                    F[xind]     += (kecm/nvtmp)*dpinx;
                    F[yind]     += (kecm/nvtmp)*dpiny;
                }
                else{
                    pinpos[NDIM*pinCellInd] = cx;
                    pinpos[NDIM*pinCellInd + 1] = cy;
                }
            }
        }
    }
}



/******************************

    F O R C E

    U P D A T E S

*******************************/


void tumor2D::resetForcesAndEnergy(){
    fill(F.begin(), F.end(), 0.0);
    fill(stress.begin(), stress.end(), 0.0);
    fill(wpress.begin(), wpress.end(), 0.0);
    U = 0.0;
    Ua=0.0;
    Ul=0.0;
    Ub=0.0;
    Utest=0.0;
}

void tumor2D::tumorShapeForces(){
    // local variables
    int ci, gi, vi, nvtmp;
    double fa, fli, flim1, fb, cx, cy, xi, yi;
    double l0im1, l0i, a0tmp, atmp;
    double dx, dy, da, dli, dlim1, dtim1, dti, dtip1;
    double lim2x, lim2y, lim1x, lim1y, lix, liy, lip1x, lip1y, li, lim1;
    double rim2x, rim2y, rim1x, rim1y, rix, riy, rip1x, rip1y, rip2x, rip2y;
    double nim1x, nim1y, nix, niy, sinim1, sini, sinip1, cosim1, cosi, cosip1;
    double ddtim1, ddti;

    // loop over vertices, add to force
    ci = 0;
    for (gi = 0; gi < NVTOT; gi++) {
        // if tumor, then no shape force
        cindices(ci,vi,gi);

        //warning
        if (nv[gi]==1) {
            continue;
        }
        // -- Area force (and get cell index ci)
        if (ci < NCELLS) {
            if (gi == szList[ci]) {
                // shape information
                nvtmp = nv[ci];
                a0tmp = a0[ci];

                // preferred segment length of last segment
                l0im1 = l0[im1[gi]];

                // compute area deviation
                atmp = area(ci);
                da = (atmp / a0tmp) - 1.0;

                // update potential energy
                U += 0.5 * ka * (da * da);
                Ua+=0.5 * ka * (da * da);
                // shape force parameters
                //warning
                //fa = ka * da * (rho0 / a0tmp);
                fa = ka * da / a0tmp;
                
                //why?
                //fb = kb * rho0;
                fb = kb;

                // compute cell center of mass
                xi = x[NDIM * gi];
                yi = x[NDIM * gi + 1];
                cx = xi;
                cy = yi;
                for (vi = 1; vi < nvtmp; vi++) {
                    // get distances between vim1 and vi
                    dx = x[NDIM * (gi + vi)] - xi;
                    dy = x[NDIM * (gi + vi) + 1] - yi;
                    if (pbc[0])
                        dx -= L[0] * round(dx / L[0]);
                    if (pbc[1])
                        dy -= L[1] * round(dy / L[1]);

                    // add to centers
                    xi += dx;
                    yi += dy;

                    cx += xi;
                    cy += yi;
                }
                cx /= nvtmp;
                cy /= nvtmp;

                // get coordinates relative to center of mass
                rix = x[NDIM * gi] - cx;
                riy = x[NDIM * gi + 1] - cy;

                // get prior adjacent vertices
                rim2x = x[NDIM * im1[im1[gi]]] - cx;
                rim2y = x[NDIM * im1[im1[gi]] + 1] - cy;
                if (pbc[0])
                    rim2x -= L[0] * round(rim2x / L[0]);
                if (pbc[1])
                    rim2y -= L[1] * round(rim2y / L[1]);

                rim1x = x[NDIM * im1[gi]] - cx;
                rim1y = x[NDIM * im1[gi] + 1] - cy;
                if (pbc[0])
                    rim1x -= L[0] * round(rim1x / L[0]);
                if (pbc[1])
                    rim1y -= L[1] * round(rim1y / L[1]);

                // get prior segment vectors
                lim2x = rim1x - rim2x;
                lim2y = rim1y - rim2y;

                lim1x = rix - rim1x;
                lim1y = riy - rim1y;

                // increment cell index
                //ci++;
            }
        }

        // preferred segment length
        l0i = l0[gi];

        // get next adjacent vertices
        rip1x = x[NDIM * ip1[gi]] - cx;
        rip1y = x[NDIM * ip1[gi] + 1] - cy;
        if (pbc[0])
            rip1x -= L[0] * round(rip1x / L[0]);
        if (pbc[1])
            rip1y -= L[1] * round(rip1y / L[1]);

        // -- Area force
        F[NDIM * gi] += 0.5 * fa * (rim1y - rip1y);
        F[NDIM * gi + 1] += 0.5 * fa * (rip1x - rim1x);

        // -- Perimeter force
        lix = rip1x - rix;
        liy = rip1y - riy;

        // segment lengths
        lim1 = sqrt(lim1x * lim1x + lim1y * lim1y);
        li = sqrt(lix * lix + liy * liy);

        // segment deviations
        dlim1 = (lim1 / l0im1) - 1.0;
        dli = (li / l0i) - 1.0;


        // segment forces
        //flim1 = kl/nv[ci] * (rho0 / l0im1);
        //fli = kl/nv[ci] * (rho0 / l0i);
        flim1 = kl/nv[ci] / l0im1;
        fli = kl/nv[ci] / l0i;

        // add to forces
        F[NDIM * gi] += (fli * dli * lix / li) - (flim1 * dlim1 * lim1x / lim1);
        F[NDIM * gi + 1] += (fli * dli * liy / li) - (flim1 * dlim1 * lim1y / lim1);

        // update potential energy
        U += 0.5 * kl/nv[ci] * (dli * dli);
        Ul +=0.5 * kl/nv[ci] * (dli * dli);
        
        // -- Bending force
        if (kb > 0 && ci > tN - 1) {
            // get ip2 for third angle
            rip2x = x[NDIM * ip1[ip1[gi]]] - cx;
            rip2y = x[NDIM * ip1[ip1[gi]] + 1] - cy;
            if (pbc[0])
                rip2x -= L[0] * round(rip2x / L[0]);
            if (pbc[1])
                rip2y -= L[1] * round(rip2y / L[1]);

            // get last segment length
            lip1x = rip2x - rip1x;
            lip1y = rip2y - rip1y;

            // get angles
            sinim1 = lim1x * lim2y - lim1y * lim2x;
            cosim1 = lim1x * lim2x + lim1y * lim2y;

            sini = lix * lim1y - liy * lim1x;
            cosi = lix * lim1x + liy * lim1y;

            sinip1 = lip1x * liy - lip1y * lix;
            cosip1 = lip1x * lix + lip1y * liy;

            // get normal vectors
            nim1x = lim1y;
            nim1y = -lim1x;

            nix = liy;
            niy = -lix;

            // get change in angles
            dtim1 = atan2(sinim1, cosim1) - t0[im1[gi]];
            dti = atan2(sini, cosi) - t0[gi];
            dtip1 = atan2(sinip1, cosip1) - t0[ip1[gi]];

            //warning: set t0
            //t0[gi] = atan2(sini, cosi);
            
            // get delta delta theta's
            ddtim1 = (dti - dtim1) / (lim1 * lim1);
            ddti = (dti - dtip1) / (li * li);

            // add to force
            F[NDIM * gi] += fb * (ddtim1 * nim1x + ddti * nix);
            F[NDIM * gi + 1] += fb * (ddtim1 * nim1y + ddti * niy);

            // update potential energy
            U += 0.5 * kb * (dti * dti);
            //Ub+=0.5 * kb * (dti * dti);
        }

        // update old coordinates
        rim2x = rim1x;
        rim1x = rix;
        rix = rip1x;

        rim2y = rim1y;
        rim1y = riy;
        riy = rip1y;

        // update old segment vectors
        lim2x = lim1x;
        lim2y = lim1y;

        lim1x = lix;
        lim1y = liy;

        // update old preferred segment length
        l0im1 = l0i;
    }
}

void tumor2D::repulsiveTumorInterfaceForces() {
    // local variables
    int ci, cj, gi, gj, vi, vj, bi, bj, pi, pj, boxid, sbtmp;
    double sij, rij, xij, dx, dy, xi, yi, ri;
    double ftmp, fx, fy;

    // sort particles
    sortNeighborLinkedList2D();

    // reset contact network
    fill(cij.begin(), cij.end(), 0);

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
                fx = kc*(1.0 - ((xi-wpos)/ri))/ri;
                F[NDIM*gi] += fx;

                // update wall stress
                wpress[0] += fx/L[1];
            }
            else if (xi > L[0] - ri){
                // update forces
                fx = -kc*(1.0 - ((L[0] - xi)/ri))/ri;
                F[NDIM*gi] += fx;

                // update wall stresses
                //wpress[0] -= fx/L[1];
            }

            // cell index of gi
            cindices(ci, vi, gi);

            // next particle in list
            pj = list[pi];

            // loop down neighbors of pi in same cell
            while (pj > 0) {
                // real index of pj
                gj = pj - 1;

                if (gj == ip1[gi] || gj == im1[gi]) {
                    pj = list[pj];
                    continue;
                }

                // cell index of gj
                cindices(cj, vj, gj);

                // contact distance
                sij = r[gi] + r[gj];

                // particle distance
                dx = x[NDIM * gj] - x[NDIM * gi];
                if (pbc[0])
                    dx -= L[0] * round(dx / L[0]);
                if (dx < sij) {
                    dy = x[NDIM * gj + 1] - x[NDIM * gi + 1];
                    if (pbc[1])
                        dy -= L[1] * round(dy / L[1]);
                    if (dy < sij) {
                        rij = sqrt(dx * dx + dy * dy);
                        if (rij < sij) {
                            // scaled distance
                            xij = rij/sij;

                            // force magnitude
                            ftmp = kc*(1 - xij)/sij;

                            // increase potential energy
                            U += 0.5*kc*pow(1.0 - xij,2.0);

                            // force elements
                            fx                     = ftmp*(dx/rij);
                            fy                     = ftmp*(dy/rij);

                            // add to forces
                            F[NDIM*gi]             -= fx;
                            F[NDIM*gi + 1]         -= fy;

                            F[NDIM*gj]             += fx;
                            F[NDIM*gj + 1]         += fy;

                            // add to virial stress
                            stress[0]             += dx*fx;
                            stress[1]             += dy*fy;
                            stress[2]             += 0.5*(dx*fy + dy*fx);

                            // add to contacts
                            if (ci > cj)
                                cij[NCELLS*cj + ci - (cj+1)*(cj+2)/2]++;
                            else if (ci < cj)
                                cij[NCELLS*ci + cj - (ci+1)*(ci+2)/2]++;
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

                    if (gj == ip1[gi] || gj == im1[gi]) {
                        pj = list[pj];
                        continue;
                    }

                    // cell index of gj
                    cindices(cj, vj, gj);

                    // contact distance
                    sij = r[gi] + r[gj];

                    // particle distance
                    dx = x[NDIM * gj] - x[NDIM * gi];
                    if (pbc[0])
                        dx -= L[0] * round(dx / L[0]);
                    if (dx < sij) {
                        dy = x[NDIM * gj + 1] - x[NDIM * gi + 1];
                        if (pbc[1])
                            dy -= L[1] * round(dy / L[1]);
                        if (dy < sij) {
                            rij = sqrt(dx * dx + dy * dy);
                            if (rij < sij) {
                                // scaled distance
                                xij = rij/sij;

                                // force magnitude
                                ftmp = kc*(1 - xij)/sij;

                                // increase potential energy
                                U += 0.5*kc*pow(1.0 - xij,2.0);

                                // force elements
                                fx                     = ftmp*(dx/rij);
                                fy                     = ftmp*(dy/rij);

                                // add to forces
                                F[NDIM*gi]             -= fx;
                                F[NDIM*gi + 1]         -= fy;

                                F[NDIM*gj]             += fx;
                                F[NDIM*gj + 1]         += fy;

                                // add to virial stress
                                stress[0]             += dx*fx;
                                stress[1]             += dy*fy;
                                stress[2]             += 0.5*(dx*fy + dy*fx);

                                // add to contacts
                                if (ci > cj)
                                    cij[NCELLS*cj + ci - (cj+1)*(cj+2)/2]++;
                                else if (ci < cj)
                                    cij[NCELLS*ci + cj - (ci+1)*(ci+2)/2]++;
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

    // normalize stress by box area, make dimensionless
    stress[0] *= (rho0 / (L[0] * L[1]));
    stress[1] *= (rho0 / (L[0] * L[1]));
    stress[2] *= (rho0 / (L[0] * L[1]));
}

void tumor2D::stickyTumorInterfaceForces(){
    // local variables
    int ci, cj, ck, gi, gj, gk, vi, vj, vk, bi, bj, pi, pj, boxid, sbtmp;
    double sij, rij, dx, dy, xi, yi, ri;
    double ftmp, fx, fy, fxtmp, fytmp;
    //miu is friction constand in the 1st method or preferred velocity in the 2nd method
    vector<int> ztt(NVTOT,0);

    // attraction shell parameters
    double shellij, cutij, xij, kint = (kc*l1)/(l2 - l1);

    // sort particles
    //sortNeighborLinkedList2D();

    // reset contact network
    fill(cij.begin(), cij.end(), 0);
    fill(F_ij.begin(),F_ij.end(),0.0);
    
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
                U += 0.5*kc*pow(1.0 - xij,2.0);
                // update wall stress
                wpress[0] += fx/L[1];
            }
            else if (xi > L[0] - ri){
                // update forces
                xij = (L[0] - xi)/ri;
                fx = -kc*(1.0 - xij)/ri;
                F[NDIM*gi] += fx;
                U += 0.5*kc*pow(1.0 - xij,2.0);

                // update wall stresses
                //wpress[0] -= fx/L[1];
            }
            
            // cell index of gi
            cindices(ci, vi, gi);
            // next particle in list
            pj = list[pi];

            // loop down neighbors of pi in same cell
            while (pj > 0) {
                // real index of pj
                gj = pj - 1;

                if (gj == ip1[gi] || gj == im1[gi]) {
                    //cout << gi << " " << gj << endl;
                    pj = list[pj];
                    continue;
                }

                // cell index of gj
                cindices(cj, vj, gj);

                // contact distance
                sij = r[gi] + r[gj];
                
                // attraction distances
                shellij = (1.0 + l2)*sij;
                cutij = (1.0 + l1)*sij;

                // particle distance
                dx = x[NDIM * gj] - x[NDIM * gi];
                if (pbc[0])
                    dx -= L[0] * round(dx / L[0]);
                if (dx < shellij) {
                    dy = x[NDIM * gj + 1] - x[NDIM * gi + 1];
                    if (pbc[1])
                        dy -= L[1] * round(dy / L[1]);
                    if (dy < shellij) {
                        rij = sqrt(dx * dx + dy * dy);
                        if (rij < shellij) {
                            // scaled distance
                            xij = rij/sij;
                            
                            // pick force based on vertex-vertex distance
                            // only tumor cells attract
                            if (rij > cutij && ci < tN && cj < tN){
                                // force scale
                                ftmp = kint*(xij - 1.0 - l2)/sij;

                                // increase potential energy
                                U += -0.5*kint*pow(1.0 + l2 - xij,2.0);
                                Utest+=-0.5*kint*pow(1.0 + l2 - xij,2.0);
                            }
                            else{
                                if ((ci < tN && cj < tN) || rij < sij){
                                    // force scale
                                    ftmp = kc*(1 - xij)/sij;

                                    // increase potential energy
                                    if(ci < tN && cj < tN){
                                        U += 0.5*kc*(pow(1.0 - xij,2.0) - l1*l2);
                                        Utest += 0.5*kc*(pow(1.0 - xij,2.0) - l1*l2);
                                    }
                                    else{
                                        U += 0.5*kc*pow(1.0 - xij,2.0);
                                    }
                                }
                                else{
                                    pj = list[pj];
                                    continue;
                                }
                            }
                            // force elements
                            fx                     = ftmp*(dx/rij);
                            fy                     = ftmp*(dy/rij);
                            
                            // add to forces
                            F[NDIM*gi]             -= fx;
                            F[NDIM*gi + 1]         -= fy;

                            F[NDIM*gj]             += fx;
                            F[NDIM*gj + 1]         += fy;

                            // add to virial stress
                            stress[0]             += dx*fx;
                            stress[1]             += dy*fy;
                            stress[2]             += 0.5*(dx*fy + dy*fx);

                            // add to contacts
                            if (ci > cj)
                                cij[NCELLS*cj + ci - (cj+1)*(cj+2)/2]++;
                            else if (ci < cj)
                                cij[NCELLS*ci + cj - (ci+1)*(ci+2)/2]++;

                            // if both tumor cells, add to contact list for surface tension
                            if (ci < tN && cj < tN){
                                ztt[gi]++;
                                ztt[gj]++;
                            }
                            
                            //if one is tumor and not attractive
                            if (xij<1) {
                                if (ci < tN){
                                    F_ij[NDIM*(ci*NCELLS + cj)]     -= fx;
                                    F_ij[NDIM*(ci*NCELLS + cj) + 1] -= fy;
                                }
                                if (cj < tN){
                                    F_ij[NDIM*(cj*NCELLS + ci)]     += fx;
                                    F_ij[NDIM*(cj*NCELLS + ci) + 1] += fy;
                                }
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

                    if (gj == ip1[gi] || gj == im1[gi]) {
                        //cout << gi << " " << gj << endl;
                        pj = list[pj];
                        continue;
                    }

                    // cell index of gj
                    cindices(cj, vj, gj);

                    // contact distance
                    sij = r[gi] + r[gj];

                    // attraction distances
                    shellij = (1.0 + l2)*sij;
                    cutij = (1.0 + l1)*sij;

                    // particle distance
                    dx = x[NDIM * gj] - x[NDIM * gi];
                    if (pbc[0])
                        dx -= L[0] * round(dx / L[0]);
                    if (dx < shellij) {
                        dy = x[NDIM * gj + 1] - x[NDIM * gi + 1];
                        if (pbc[1])
                            dy -= L[1] * round(dy / L[1]);
                        if (dy < shellij) {
                            rij = sqrt(dx * dx + dy * dy);
                            if (rij < shellij) {
                                // scaled distance
                                xij = rij/sij;

                                // pick force based on vertex-vertex distance
                                // only tumor cells attract
                                if (rij > cutij && ci < tN && cj < tN){
                                    // force scale
                                    ftmp = kint*(xij - 1.0 - l2)/sij;

                                    // increase potential energy
                                    U += -0.5*kint*pow(1.0 + l2 - xij,2.0);
                                    Utest+=-0.5*kint*pow(1.0 + l2 - xij,2.0);
                                }
                                else{
                                    if ((ci < tN && cj < tN) || rij < sij){
                                        // force scale
                                        ftmp = kc*(1 - xij)/sij;

                                        // increase potential energy
                                        if(ci < tN && cj < tN){
                                            U += 0.5*kc*(pow(1.0 - xij,2.0) - l1*l2);
                                            Utest += 0.5*kc*(pow(1.0 - xij,2.0) - l1*l2);
                                        }
                                        else{
                                            U += 0.5*kc*pow(1.0 - xij,2.0);
                                        }
                                        
                                    }
                                    else{
                                        pj = list[pj];
                                        continue;
                                    }
                                }
                                // force elements
                                fx                     = ftmp*(dx/rij);
                                fy                     = ftmp*(dy/rij);

                                
                                // add to forces
                                F[NDIM*gi]             -= fx;
                                F[NDIM*gi + 1]         -= fy;

                                F[NDIM*gj]             += fx;
                                F[NDIM*gj + 1]         += fy;

                                // add to virial stress
                                stress[0]             += dx*fx;
                                stress[1]             += dy*fy;
                                stress[2]             += 0.5*(dx*fy + dy*fx);
                                
                                // add to contacts
                                if (ci > cj)
                                    cij[NCELLS*cj + ci - (cj+1)*(cj+2)/2]++;
                                else if (ci < cj)
                                    cij[NCELLS*ci + cj - (ci+1)*(ci+2)/2]++;

                                // if both tumor cells, add to contact list for surface tension
                                if (ci < tN && cj < tN){
                                    ztt[gi]++;
                                    ztt[gj]++;
                                }
                                
                                
                                //if one is tumor and not attractive
                                if (xij<1) {
                                    if (ci < tN){
                                        F_ij[NDIM*(ci*NCELLS + cj)]     -= fx;
                                        F_ij[NDIM*(ci*NCELLS + cj) + 1] -= fy;
                                    }
                                    if (cj < tN){
                                        F_ij[NDIM*(cj*NCELLS + ci)]     += fx;
                                        F_ij[NDIM*(cj*NCELLS + ci) + 1] += fy;
                                    }
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
    
    //update contactTime
    for (ci=0; ci<tN; ci++){
        for (cj=tN; cj<NCELLS; cj++){
            if (F_ij[NDIM*(ci*NCELLS + cj)] == 0.0) {
                contactTime[ci * (NCELLS - tN) + cj - tN] = 0;
            }
            else {
                contactTime[ci * (NCELLS - tN) + cj - tN] ++;
            }
        }
    }

    
    // normalize stress by box area, make dimensionless
    stress[0] *= (rho0 / (L[0] * L[1]));
    stress[1] *= (rho0 / (L[0] * L[1]));
    stress[2] *= (rho0 / (L[0] * L[1]));
    // update l0 based on ztt
    /*
    for (gi=0; gi<NVTOT; gi++){
        if (ztt[gi] == 0 && ztt[ip1[gi]] == 0)
            l0[gi] = l0_init[gi]*(1.0 - (gamtt/kl));
        else
            l0[gi] = l0_init[gi];
    }
     */
}

void tumor2D::repulsiveTumorInterfaceForceUpdate() {
    resetForcesAndEnergy();
    repulsiveTumorInterfaceForces();
    tumorShapeForces();
}

void tumor2D::stickyTumorInterfaceForceUpdate() {
    resetForcesAndEnergy();
    //crawlerUpdate();
    stickyTumorInterfaceForces();
    tumorShapeForces();
    adipocyteECMAdhesionForces();
}





/******************************

    T U M O R

    F I R E

*******************************/


void tumor2D::tumorFIRE(tumor2DMemFn forceCall, double Ftol, double dt0) {
    // local variables
    int i;

    // check to see if cell linked-list has been initialized
    if (NBX == -1) {
        cerr << "    ** ERROR: In dpm::fire, NBX = -1, so cell linked-list has not yet been initialized. Ending here.\n";
        exit(1);
    }

    // FIRE variables
    double P, fnorm, fcheck, vnorm, alpha, dtmax, dtmin, frec;
    int npPos, npNeg, fireit;

    // set dt based on geometric parameters
    setdt(dt0);

    // Initialize FIRE variables
    P = 0;
    fnorm = 0;
    vnorm = 0;
    alpha = alpha0;

    dtmax = 10.0 * dt;
    //1e-2
    dtmin = 0;

    npPos = 0;
    npNeg = 0;

    fireit = 0;
    fcheck = 10 * Ftol;
    frec = 0.0;
    
    // reset forces and velocities
    resetForcesAndEnergy();
    fill(v.begin(), v.end(), 0.0);

    // relax forces using FIRE
    while ((fcheck > Ftol || fireit < NDELAY) && fireit < itmax) {
        // compute P
        P = 0.0;
        for (i = 0; i < vertDOF; i++)
            P += v[i] * F[i];

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
            cout << "    ** U         = " << U << endl;
            cout << "    ** dt         = " << dt << endl;
            cout << "    ** P         = " << P << endl;
            cout << "    ** alpha     = " << alpha << endl;
            cout << "    ** npPos     = " << npPos << endl;
            cout << "   ** npNeg     = " << npNeg << endl;
            cout << "   ** phi      = " << vertexPreferredPackingFraction2D() << endl;
            
            if (fcheck > 0.99 * frec && fcheck < 1.01 * frec) {
                break;
            }
            else{
                frec = fcheck;
            }
        }

        //cout << accumulate(F.begin(), F.end(), 1) << endl;
        //cout << vnorm << endl;
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
                x[i] -= 0.5 * dt * v[i];

                // reset vertex velocities
                v[i] = 0.0;
            }

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

        // compute fnorm, vnorm and P
        fnorm = 0.0;
        vnorm = 0.0;
        for (i = 0; i < vertDOF; i++) {
            fnorm += F[i] * F[i];
            vnorm += v[i] * v[i];
        }
        fnorm = sqrt(fnorm);
        vnorm = sqrt(vnorm);

        // update velocities (s.d. vs inertial dynamics) only if forces are acting
        if (fnorm > 0) {
            for (i = 0; i < vertDOF; i++)
                v[i] = (1 - alpha) * v[i] + alpha * (F[i] / fnorm) * vnorm;
        }

        // VV POSITION UPDATE
        for (i = 0; i < vertDOF; i++) {
            // update position
            x[i] += dt * v[i];

            // recenter in box
            if (x[i] > L[i % NDIM] && pbc[i % NDIM])
                x[i] -= L[i % NDIM];
            else if (x[i] < 0.0 && pbc[i % NDIM])
                x[i] += L[i % NDIM];
        }

        
        // update forces (function passed as argument)
        CALL_MEMBER_FN(*this, forceCall)();

        // VV VELOCITY UPDATE #2
        for (i = 0; i < vertDOF; i++){
            v[i] += 0.5 * F[i] * dt;
        }

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
        // update iterator
        fireit++;

    }
    // check if FIRE converged
    if (fireit == itmax) {
        cout << "    ** FIRE minimization did not converge, fireit = " << fireit << ", itmax = " << itmax << "; ending." << endl;
        //warning
        exit(1);
    }
    else {
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
    }
}



/******************************

    P R O T O C O L S

*******************************/


void tumor2D::setupCheck(){
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


    // check initialization of shape parameters
    if (l0.size() != NVTOT){
        cerr << "** ERROR: in setupCheck, l0.size = " << l0.size() << ", which is != NVTOT = " << NVTOT << ". Ending. " << endl;
        exit(1);
    }
    if (l0_init.size() != NVTOT){
        cerr << "** ERROR: in setupCheck, l0_init.size = " << l0_init.size() << ", which is != NVTOT = " << NVTOT << ". Ending. " << endl;
        exit(1);
    }
    if (t0.size() != NVTOT){
        cerr << "** ERROR: in setupCheck, t0.size = " << t0.size() << ", which is != NVTOT = " << NVTOT << ". Ending. " << endl;
        exit(1);
    }
    if (a0.size() != NCELLS){
        cerr << "** ERROR: in setupCheck, a0.size = " << a0.size() << ", which is != NCELLS = " << NCELLS << ". Ending. " << endl;
        exit(1);
    }
}


// compression with boundary forces
void tumor2D::tumorCompression(double Ftol, double Ptol, double dt0, double dphi0){
    int ci;
    double meanArea=0.0;
    // check correct setup
    setupCheck();

    // local variables
    int k = 0, nr;
    double pcheck=0.0, phi0, scaleFactor = 1.0;

    // initialize preferred packing fraction
    phi0 = vertexPreferredPackingFraction2D();

    // loop until pcheck > Ptol is found
    //warning
    while (pcheck < Ptol && k < itmax) {
        // relax configuration (pass repsulive force update member function)
        // scale particle sizes
        // update packing fraction
        phi0 = vertexPreferredPackingFraction2D();

        scaleFactor = sqrt((phi0 + dphi0) / phi0);
        scaleParticleSizes2D(scaleFactor);
        // rescale so that mean(a0(tN+1,NCELLS))=1
        //scaleFactor=tumorRescale(dt0);
        //Ftol = Ftol * scaleFactor;
        tumorFIRE(&tumor2D::repulsiveTumorInterfaceForceUpdate, Ftol, dt0);

        // update pressure
        //pcheck = 0.5 * (stress[0] + stress[1]);
        rho0 = 0;
        for(ci=tN;ci<NCELLS;ci++){
            rho0 +=a0[ci];
        }
        rho0 = rho0/(NCELLS-tN);
        //warning
        rho0 = pow((r.at(1)+r.at(800))/2/0.070523697943470,2.0);
        pcheck = wpress[0] * rho0;
        
        // remove rattlers
        nr = removeRattlers();

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
        cout << "    * phi             = " << vertexPackingFraction2D() << endl;
        cout << "    * scaleFactor     = " << scaleFactor << endl;
        cout << "    * pcheck         = " << pcheck << endl;
        cout << "    * U              = " << U << endl;
        cout << "    * Nvv              = " << vvContacts() << endl;
        cout << "    * Ncc             = " << ccContacts() << endl;
        cout << endl << endl;
        
    
        printTumorInterface(0.0);
        // update iterate
        k++;
    }
    scaleFactor=tumorRescale(dt0);
    fill(v.begin(), v.end(), 0.0);
    for (ci=0; ci<tN; ci++){
        psi[ci] = 2.0*PI*drand48();
    }
}

//rescale so that mean(a0(tN+1,NCELLS))=1
double tumor2D::tumorRescale(double dt0){
    int ci, gi, vi;
    double meanArea=0.0, scaleFactor=1.0;
    //calculate rho0
    for(ci=tN;ci<NCELLS;ci++){
        meanArea +=a0[ci];
    }
    meanArea = meanArea/(NCELLS-tN);
    
    scaleFactor = sqrt(meanArea);
    //make r0=0.070523697943470;
    scaleFactor = (r.at(1)+r.at(800))/2/0.070523697943470;
    
    L[0] = L[0]/scaleFactor;
    L[1] = L[1]/scaleFactor;
    
    for (ci = 0; ci < NCELLS; ci++) {
        a0[ci] = a0[ci]/scaleFactor/scaleFactor;
        
        gi = gindex(ci, 0);
        x.at(NDIM * gi) = x.at(NDIM * gi)/scaleFactor;
        x.at(NDIM * gi + 1) = x.at(NDIM * gi + 1)/scaleFactor;
        r.at(gi) = r.at(gi)/scaleFactor;
        l0.at(gi) = l0.at(gi)/scaleFactor;
        for (vi = 1; vi < nv.at(ci); vi++) {
            gi++;
            x.at(NDIM * gi) = x.at(NDIM * gi)/scaleFactor;
            x.at(NDIM * gi + 1) = x.at(NDIM * gi + 1)/scaleFactor;
            r.at(gi) = r.at(gi)/scaleFactor;
            l0.at(gi) = l0.at(gi)/scaleFactor;
        }
    }
    
    setdt(dt0);
    
    return scaleFactor;
}

// invasion protocol
void tumor2D::invasion(tumor2DMemFn forceCall, double dDr, double dPsi, double Drmin, int NT, int NPRINTSKIP){
    // check correct setup
    setupCheck();

    // local variables
    int k, i, ci, cj;
    double t = 0.0, zta, Drtmp;

    // attach pins
    updateECMAttachments(1);

    //initialize psi
    for (ci=0; ci<tN; ci++){
        psi[ci] = 2.0*PI*drand48();
    }
    // loop over time, have active brownian crawlers invade adipocytes
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
        
        //adipocyte shrink response to pressure
        adipocyteShrink();
        
        // update forces
        CALL_MEMBER_FN(*this, forceCall)();
        
        //update psi
        psiECM();
        
        // update psi based on persistence
        psiDiffusion();
        
        // update active brownian crawler
        //crawlerUpdate();
                
        // update positions (EULER UPDATE, OVERDAMPED)
        for (i=0; i<vertDOF; i++)
            x[i] += dt * F[i];

        // update time
        t += dt;

        // print message console, print position to file
        if (k % NPRINTSKIP == 0){
            cout << endl << endl;
            cout << "===========================================" << endl;
            cout << "            invading tumor cells             " << endl;
            cout << "===========================================" << endl;
            cout << endl;
            cout << "    ** k             = " << k << endl;
            cout << "    ** p             = " << wpress[0] << endl;
            cout << "    ** phi             = " << vertexPackingFraction2D() << endl;

            // print vertex positions to check placement
            cout << "\t** PRINTING POSITIONS TO FILE... " << endl;
            printTumorInterface(t);
        }
    }
}

// invasion at constant pressure
void tumor2D::invasionConstP(tumor2DMemFn forceCall, double M, double P0, double g0, double dDr, double dPsi, double Drmin, int NT, int NPRINTSKIP){
    // check correct setup
    setupCheck();
    
    // local variables
    int k, i, ci, cj, gi, d, upb = 20000;
    double t = 0.0, zta, Drtmp, Lold, Lnew, gam;
    double subBoxLength = 2.0;
    double x_max=0;
    int press_teller,press_it=0;
    double press_ave =0.0;
    double dw = 0.0;
    double B = 0.0;
    double M_wall = M*tN;
    double V_wall = 0.0;
    double K_t=0.0;
    double K=0.0;
    
    vector<int> tN_list;
    for (ci=0; ci<tN; ci++){
        tN_list.push_back(ci);
    }
    random_shuffle (tN_list.begin(), tN_list.end());

    // attach pins
    updateECMAttachments(1);

    //initialize psi warning
    for (ci=0; ci<tN; ci++){
        psi[ci] = 2.0*PI*tN_list[ci]/tN;
    }
    //initialize velocity
    for (ci=0; ci<tN; ci++){
        v[NDIM*ci] = v0*cos(psi[ci]);
        v[NDIM*ci + 1] = v0*sin(psi[ci]);
    }
    
    
    //fill(v.begin(), v.end(), 0.0);
    fill(t0.begin(), t0.end(), 0.0);
    
    Lold = L[0];
    Lnew = Lold;
    //allocation enough memory
    x.reserve(upb);
    v.reserve(upb);
    F.reserve(upb);
    nv.reserve(upb);
    r.reserve(upb);
    a0.reserve(upb);
    l0.reserve(upb);
    t0.reserve(upb);
    l0_init.reserve(upb);
    Dr.reserve(upb);
    psi.reserve(upb);
    cij.reserve(40000);
    F_ij.reserve(40000);
    contactTime.reserve(40000);
    list.reserve(upb);
    szList.reserve(upb);
    pinattach.reserve(upb);
    ifbroken.reserve(upb);
    
    reNeighborLinkedList2D(subBoxLength);
    neighborLinkedList2D();
    // initial pressure
    CALL_MEMBER_FN(*this, forceCall)();
    
    // RELAXATION: update tumor cell positions (EULER UPDATE, OVERDAMPED)
    press_teller = 1;
    while (press_teller ==0) {
        press_it += 1;
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
        for (i=0; i<tN*NDIM; i++){
            x[i] += dt * (v[i] +dt/2.0/M * (F[i] -B*v[i]));
            v[i] += dt/2.0/M * (F[i] -B*v[i]*  (1.0+1.0/(1.0+B/2.0*dt)));
        }
        wpos += dt*(V_wall + dt/2.0/M_wall*((P0 - wpress[0])*L[1]-B*V_wall));
        V_wall += dt/2.0/M_wall * ((P0 - wpress[0])*L[1]-B*V_wall*(1.0+1.0/(1.0+B/2.0*dt)));
        
        
        // update psi before update force
        // update psi based on persistence warning
        psiDiffusion();
        //update psi warning
        //psiECM();
        // update forces
        reNeighborLinkedList2D(subBoxLength);
        neighborLinkedList2D();
        resetForcesAndEnergy();
        //crawlerUpdate();
        stickyTumorInterfaceForces();


        // update velocity 2nd term (Velocity Verlet, OVERDAMPED)
        for (i=0; i<tN*NDIM; i++)
            v[i] += dt/2.0/M * F[i]/(1.0+B/2.0*dt);

        V_wall += dt/2.0/M_wall * (P0 - wpress[0])*L[1] / (1.0+B/2.0*dt);
      
        /*******************************************************************************************************************************/
        
        if (press_it % 10000 == 0) {
            //kinetic energy
            K=0;
            K_t=0;
            for (int i = 0; i < NDIM*tN; i++){
                K += v[i] * v[i];
                K_t += v[i] * v[i];
            }
            for (int i = NDIM*tN; i < vertDOF; i++)
                K += v[i] * v[i];
            K *= 0.5;
            K_t *= 0.5;
            
            cout << "===========================================" << endl;
            cout << "            Relaxation                     " << endl;
            cout << "===========================================" << endl;
            cout << endl;
            cout << "    ** press_it              = " << press_it << endl;
            cout << "    ** press_ave             = " << press_ave << endl;
            cout << "    ** P0                    = " << P0 << endl;
            cout << "    ** wpos                  = " << wpos << endl;
            cout << "   ** E              = " << U+K_t - wpos*P0*L[1] << endl;
            cout << "   ** U              = " << U << endl;
            cout << "   ** Kinetic        = " << K << endl;
            cout << "   ** Kinetic_tumor  = " << K_t << endl;
            cout << "   ** Kinetic_wall   = " << M_wall*V_wall*V_wall/2 << endl;
            
        }
        
        if(press_it>20000){
            press_ave += wpress[0];
            if (press_it % 10000 == 0) {
                press_ave = press_ave / 10000;
                if (press_ave < 1.01* P0 && press_ave > 0.99* P0) {
                    press_teller = 1;
                } else {
                    press_ave = 0;
                }
            }
        }
    }

    fill(contactTime.begin(),contactTime.end(),0);
    // update forces
    // loop over time, have active brownian crawlers invade adipocytes
    
    //printTumorInterface(t);
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
            x[i] += dt * (v[i] +dt/2.0/M * (F[i] -B*v[i]));
            v[i] += dt/2.0/M * (F[i] -B*v[i]*  (1.0+1.0/(1.0+B/2.0*dt)));
        }
        wpos += dt*(V_wall + dt/2.0/M_wall*((P0 - wpress[0])*L[1]-B*V_wall));
        V_wall += dt/2.0/M_wall * ((P0 - wpress[0])*L[1]-B*V_wall*(1.0+1.0/(1.0+B/2.0*dt)));
        // update psi before update force
        psiDiffusion();
        //psiECM();
        
        // sort particles
        reNeighborLinkedList2D(subBoxLength);
        neighborLinkedList2D();
        // update forces
        CALL_MEMBER_FN(*this, forceCall)();
        /*//Euler update
        for (i=0; i<vertDOF; i++){
            v[i] += dt * F[i]/M;
            x[i] += dt * v[i];
        }
        V_wall += dt * (P0 - wpress[0])*L[1] / M_wall;
        wpos += dt * V_wall;
         */
        // update velocity 2nd term (Velocity Verlet, OVERDAMPED)
        for (i=0; i<vertDOF; i++)
            v[i] += dt/2.0/M * F[i]/(1.0+B/2.0*dt);

        V_wall += dt/2.0/M_wall * (P0 - wpress[0])*L[1] / (1.0+B/2.0*dt);

        
        //kinetic energy
        K=0;
        K_t=0;
        for (int i = 0; i < NDIM*tN; i++){
            K += v[i] * v[i];
            K_t += v[i] * v[i];
        }
        for (int i = NDIM*tN; i < vertDOF; i++)
            K += v[i] * v[i];
        K *= 0.5;
        K_t *= 0.5;
        
        
        if(k<100000){
            for (int i = 0; i < vertDOF; i++){
                v[i] = v[i] * sqrt(NVTOT*M*v0*v0/2/K);
            }
            //kinetic energy
            K=0;
            K_t=0;
            for (int i = 0; i < NDIM*tN; i++){
                K += v[i] * v[i];
                K_t += v[i] * v[i];
            }
            for (int i = NDIM*tN; i < vertDOF; i++)
                K += v[i] * v[i];
            K *= 0.5;
            K_t *= 0.5;
        }
        
        
        // update time
        t += dt;
        /********************************************************************************************************************************/
        
        //tumor cell growth and divide warning
        //tumorGrowth(g0);
        //tumorDivide(g0);
        
        // print message console, print position to file
        if ((k+1) % NPRINTSKIP == 0 || k==0){
            
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
            cout << "   ** N_t            = " << tN << endl;
            cout << "    ** P wall        = " << wpress[0] << endl;
            cout << "    ** Lx            = " << L[0]-wpos << endl;
            cout << "   ** front          = " << x_max << endl;
            cout << "   ** E              = " << U + K + M_wall*V_wall*V_wall/2 - wpos*P0*L[1]  << endl;
            cout << "   ** U              = " << U << endl;
            cout << "   ** Ua             = " << Ua << endl;
            cout << "   ** Ul             = " << Ul << endl;
            cout << "   ** Kinetic        = " << K<< endl;
            cout << "   ** Kinetic_tumor  = " << K_t << endl;
            cout << "   ** Kinetic_wall   = " << M_wall*V_wall*V_wall/2 << endl;
            cout << "   ** potential_wall = " << -wpos*P0*L[1] << endl;
            cout << "   ** U_tumor        = " << Utest << endl;
            //cout << "   ** rij              = " << sqrt((x[0]-x[2])*(x[0]-x[2])+(x[1]-x[3])*(x[1]-x[3])) << endl;
            //cout << "   ** Fij              = " << sqrt(F[0]*F[0]+F[1]*F[1]) << endl;
            // print vertex positions to check placement
            cout << "\t** PRINTING POSITIONS TO FILE... " << endl;
            
            if ((k+1) % (NPRINTSKIP*10) == 0) {
                printTumorInterface(t);
            }
        }
    }
}


// just get them bois crawlin'
// always print unbox coordinates for better MSD computation
void tumor2D::crawling(tumor2DMemFn forceCall, tumor2DMemFn psiCall, int NT, int NPRINTSKIP){
    // check correct setup
    setupCheck();

    // local variables
    int k, i, ci, cj;
    double t = 0.0;

    // loop over time, have cells crawl around
    for (k=0; k<NT; k++){
        // // pbcs
        // for (i=0; i<vertDOF; i++){
        //     if (x[i] > L[i % NDIM] && pbc[i % NDIM])
        //         x[i] -= L[i % NDIM];
        //     else if (x[i] < 0.0 && pbc[i % NDIM])
        //         x[i] += L[i % NDIM];
        // }

        // update forces
        CALL_MEMBER_FN(*this, forceCall)();

        // update crawling forces
        crawlerUpdate();

        // update positions (EULER UPDATE, OVERDAMPED)
        for (i=0; i<vertDOF; i++)
            x[i] += dt * F[i];

        // update psi based on persistence
        CALL_MEMBER_FN(*this, psiCall)();

        // update time
        t += dt;

        // print message console, print position to file
        if (k % NPRINTSKIP == 0){
            cout << endl << endl;
            cout << "===========================================" << endl;
            cout << "            active tumor cells                 " << endl;
            cout << "===========================================" << endl;
            cout << endl;
            cout << "    ** k             = " << k << endl;
            cout << "    ** p             = " << wpress[0] << endl;
            cout << "    ** phi             = " << vertexPackingFraction2D() << endl;

            // print vertex positions to check placement
            cout << "\t** PRINTING POSITIONS TO FILE... " << endl;
            printTumorCells(t);
        }
    }
}





/******************************

    P R I N T I N G

    F U N C T I O N S

*******************************/


void tumor2D::printTumorInterface(double t){
    // local variables
    int ci, cj, vi, gi, ctmp, zc, zv;
    double xi, yi, dx, dy, Lx, Ly;

    // check if pos object is open
    if (!posout.is_open()) {
        cerr << "** ERROR: in printConfiguration2D, posout is not open, but function call will try to use. Ending here." << endl;
        exit(1);
    }
    else
        cout << "** In printConfiguration2D, printing particle positions to file..." << endl;

    // save box sizes
    Lx = L.at(0);
    Ly = L.at(1);

    // print information starting information
    posout << setw(w) << left << "NEWFR"
           << " " << endl;
    posout << setw(w) << left << "NUMCL" << setw(w) << left << NCELLS << setw(w) << left << tN << endl;
    posout << setw(w) << left << "TSTEP" << setw(wnum) << setprecision(pnum) << left << dt << endl;
    posout << setw(w) << left << "TCURR" << setw(wnum) << setprecision(pnum) << left << t << endl;
    posout << setw(w) << left << "PACKF" << setw(wnum) << setprecision(pnum) << left << vertexPackingFraction2D() << endl;

    // print box sizes
    posout << setw(w) << left << "BOXSZ";
    posout << setw(wnum) << setprecision(pnum) << left << Lx;
    posout << setw(wnum) << setprecision(pnum) << left << Ly;
    posout << setw(wnum) << setprecision(pnum) << left << wpos;
    posout << endl;

    // print stress info
    posout << setw(w) << left << "STRSS";
    posout << setw(wnum) << setprecision(pnum) << left << stress.at(0);
    posout << setw(wnum) << setprecision(pnum) << left << stress.at(1);
    posout << setw(wnum) << setprecision(pnum) << left << stress.at(2);
    posout << endl;

    // print wall stress info
    posout << setw(w) << left << "WPRSS";
    posout << setw(wnum) << setprecision(pnum) << left << wpress.at(0);
    posout << setw(wnum) << setprecision(pnum) << left << wpress.at(1);
    posout << endl;

    // print coordinate for rest of the cells
    for (ci = 0; ci < NCELLS; ci++) {
        // get cell contact data
        zc = 0;
        zv = 0;
        for (cj = 0; cj < NCELLS; cj++) {
            if (ci != cj) {
                // contact info from entry ci, cj
                if (ci < cj)
                    ctmp = cij[NCELLS * ci + cj - (ci + 1) * (ci + 2) / 2];
                else
                    ctmp = cij[NCELLS * cj + ci - (cj + 1) * (cj + 2) / 2];

                // add to contact information
                zv += ctmp;
                if (ctmp > 0)
                    zc++;
            }
        }

        // cell information
        posout << setw(w) << left << "CINFO";
        posout << setw(w) << left << nv.at(ci);
        posout << setw(w) << left << zc;
        posout << setw(w) << left << zv;
        posout << setw(wnum) << left << a0.at(ci);
        posout << setw(wnum) << left << area(ci);
        posout << setw(wnum) << left << perimeter(ci);
        if (ci < tN){
            posout << setw(wnum) << left << psi.at(ci);
            posout << setw(wnum) << left << Dr.at(ci);
        }
        else{
            posout << setw(wnum) << left << pinpos[NDIM*(ci-tN)];
            posout << setw(wnum) << left << pinpos[NDIM*(ci-tN) + 1];
            posout << setw(w) << left << pinattach[ci-tN];
        }
        posout << endl;

        // get initial vertex positions
        gi = gindex(ci, 0);
        xi = x.at(NDIM * gi);
        yi = x.at(NDIM * gi + 1);

        // place back in box center
        if (pbc[0])
            xi = fmod(xi, Lx);
        if (pbc[1])
            yi = fmod(yi, Ly);

        posout << setw(w) << left << "VINFO";
        posout << setw(w) << left << ci;
        posout << setw(w) << left << 0;

        // output initial vertex information
        posout << setw(wnum) << setprecision(pnum) << right << xi;
        posout << setw(wnum) << setprecision(pnum) << right << yi;
        posout << setw(wnum) << setprecision(pnum) << right << r.at(gi);
        posout << setw(wnum) << setprecision(pnum) << right << l0.at(gi);
        posout << setw(wnum) << setprecision(pnum) << right << t0.at(gi);
        posout << setw(wnum) << setprecision(pnum) << right << v.at(NDIM * gi);
        posout << setw(wnum) << setprecision(pnum) << right << v.at(NDIM * gi + 1);
        posout << endl;

        // vertex information for next vertices
        for (vi = 1; vi < nv.at(ci); vi++) {
            // get global vertex index for next vertex
            gi++;

            // get next vertex positions
            dx = x.at(NDIM * gi) - xi;
            if (pbc[0])
                dx -= Lx * round(dx / Lx);
            xi += dx;

            dy = x.at(NDIM * gi + 1) - yi;
            if (pbc[1])
                dy -= Ly * round(dy / Ly);
            yi += dy;

            // Print indexing information
            posout << setw(w) << left << "VINFO";
            posout << setw(w) << left << ci;
            posout << setw(w) << left << vi;

            // output vertex information
            posout << setw(wnum) << setprecision(pnum) << right << xi;
            posout << setw(wnum) << setprecision(pnum) << right << yi;
            posout << setw(wnum) << setprecision(pnum) << right << r.at(gi);
            posout << setw(wnum) << setprecision(pnum) << right << l0.at(gi);
            posout << setw(wnum) << setprecision(pnum) << right << t0.at(gi);
            posout << setw(wnum) << setprecision(pnum) << right << v.at(NDIM * gi);
            posout << setw(wnum) << setprecision(pnum) << right << v.at(NDIM * gi + 1);
            posout << endl;
        }
    }

    // print end frame
    posout << setw(w) << left << "ENDFR" << " " << endl;
}

void tumor2D::printTumorCells(double t){
    // local variables
    int ci, cj, vi, gi, ctmp, zc, zv;
    double xi, yi, dx, dy, Lx, Ly;

    // check if pos object is open
    if (!posout.is_open()) {
        cerr << "** ERROR: in printConfiguration2D, posout is not open, but function call will try to use. Ending here." << endl;
        exit(1);
    }
    else
        cout << "** In printConfiguration2D, printing particle positions to file..." << endl;

    // save box sizes
    Lx = L.at(0);
    Ly = L.at(1);

    // print information starting information
    posout << setw(w) << left << "NEWFR"
           << " " << endl;
    posout << setw(w) << left << "NUMCL" << setw(w) << left << NCELLS << setw(w) << left << tN << endl;
    posout << setw(w) << left << "TSTEP" << setw(wnum) << setprecision(pnum) << left << dt << endl;
    posout << setw(w) << left << "TCURR" << setw(wnum) << setprecision(pnum) << left << t << endl;
    posout << setw(w) << left << "PACKF" << setw(wnum) << setprecision(pnum) << left << vertexPackingFraction2D() << endl;

    // print box sizes
    posout << setw(w) << left << "BOXSZ";
    posout << setw(wnum) << setprecision(pnum) << left << Lx;
    posout << setw(wnum) << setprecision(pnum) << left << Ly;
    posout << setw(wnum) << setprecision(pnum) << left << wpos;
    posout << endl;

    // print stress info
    posout << setw(w) << left << "STRSS";
    posout << setw(wnum) << setprecision(pnum) << left << stress.at(0);
    posout << setw(wnum) << setprecision(pnum) << left << stress.at(1);
    posout << setw(wnum) << setprecision(pnum) << left << stress.at(2);
    posout << endl;

    // print coordinate for rest of the cells
    for (ci = 0; ci < NCELLS; ci++) {
        // get cell contact data
        zc = 0;
        zv = 0;
        for (cj = 0; cj < NCELLS; cj++) {
            if (ci != cj) {
                // contact info from entry ci, cj
                if (ci < cj)
                    ctmp = cij[NCELLS * ci + cj - (ci + 1) * (ci + 2) / 2];
                else
                    ctmp = cij[NCELLS * cj + ci - (cj + 1) * (cj + 2) / 2];

                // add to contact information
                zv += ctmp;
                if (ctmp > 0)
                    zc++;
            }
        }

        // cell information
        posout << setw(w) << left << "CINFO";
        posout << setw(w) << left << nv.at(ci);
        posout << setw(w) << left << zc;
        posout << setw(w) << left << zv;
        posout << setw(wnum) << left << a0.at(ci);
        posout << setw(wnum) << left << area(ci);
        posout << setw(wnum) << left << perimeter(ci);
        posout << setw(wnum) << left << psi.at(ci);
        posout << setw(wnum) << left << Dr.at(ci);
        posout << endl;

        // get initial vertex positions
        gi = gindex(ci, 0);
        xi = x.at(NDIM * gi);
        yi = x.at(NDIM * gi + 1);

        // place back in box center
        if (pbc[0])
            xi = fmod(xi, Lx);
        if (pbc[1])
            yi = fmod(yi, Ly);

        posout << setw(w) << left << "VINFO";
        posout << setw(w) << left << ci;
        posout << setw(w) << left << 0;

        // output initial vertex information
        posout << setw(wnum) << setprecision(pnum) << right << xi;
        posout << setw(wnum) << setprecision(pnum) << right << yi;
        posout << setw(wnum) << setprecision(pnum) << right << r.at(gi);
        posout << setw(wnum) << setprecision(pnum) << right << l0.at(gi);
        posout << setw(wnum) << setprecision(pnum) << right << t0.at(gi);
        posout << endl;

        // vertex information for next vertices
        for (vi = 1; vi < nv.at(ci); vi++) {
            // get global vertex index for next vertex
            gi++;

            // get next vertex positions
            dx = x.at(NDIM * gi) - xi;
            if (pbc[0])
                dx -= Lx * round(dx / Lx);
            xi += dx;

            dy = x.at(NDIM * gi + 1) - yi;
            if (pbc[1])
                dy -= Ly * round(dy / Ly);
            yi += dy;

            // Print indexing information
            posout << setw(w) << left << "VINFO";
            posout << setw(w) << left << ci;
            posout << setw(w) << left << vi;

            // output vertex information
            posout << setw(wnum) << setprecision(pnum) << right << xi;
            posout << setw(wnum) << setprecision(pnum) << right << yi;
            posout << setw(wnum) << setprecision(pnum) << right << r.at(gi);
            posout << setw(wnum) << setprecision(pnum) << right << l0.at(gi);
            posout << setw(wnum) << setprecision(pnum) << right << t0.at(gi);
            posout << endl;
        }
    }

    // print end frame
    posout << setw(w) << left << "ENDFR" << " " << endl;
}







