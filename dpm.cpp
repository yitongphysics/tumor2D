#include "dpm.h"
#include <functional>

// namespace
using namespace std;

/******************************

	C O N S T R U C T O R S  & 

		D E S T R U C T O R

*******************************/

// Constructor with only one input (dim)
dpm::dpm(int ndim) {
	// local variables
	int d, i;

	// print to console
	cout << "** Instantiating dpm object, ndim = " << ndim << endl;

	// main variables
	NDIM = ndim;
	NNN = 13;

	// set scalars to default values
	dt = 0.0;

	kv = 0.0;
	ka = 0.0;
	kb = 0.0;
	kc = 0.0;

	l1 = 0.0;
	l2 = 0.0;

	// default boundary variables
	L.resize(NDIM);
	pbc.resize(NDIM);
	for (d = 0; d < NDIM; d++) {
		L[d] = 1.0;
		pbc[d] = 1;
	}

	// initialize nearest neighbor info
	NBX = -1;
}

// Main constructor
dpm::dpm(int n, int ndim, int seed) {
	// local variables
	int d, i;

	// print to console
	cout << "** Instantiating dpm object, NCELLS = " << n << ",  ndim = " << ndim << ", seed = " << seed << " ..." << endl;

	// main variables
	NCELLS = n;
	NDIM = ndim;
    NNN = 13;

    //polyhedron
    std::vector<std::vector<int>>        edgelist(ENUM, vector<int>(6, 0));
    std::vector<std::vector<int>>        f_unit(FNUM, vector<int>(3, 0));
    std::vector<std::vector<int>>        nnlist(NN_NUM, vector<int>(2, 0));
    std::vector<double>             A0_unit(FNUM, 0.0), A0(FNUM, 0.0);
    std::vector<double>             xyz_unit(3 * SNUM, 0.0);
    
    D0_inter_unit = 0.0;
	// set scalars to default values
	dt = 0.0;

	kv = 0.0;
	ka = 0.0;
	kb = 0.0;
	kc = 0.0;

	l1 = 0.0;
	l2 = 0.0;

	// default boundary variables
	L.resize(NDIM);
	pbc.resize(NDIM);
	for (d = 0; d < NDIM; d++) {
		L[d] = 1.0;
		pbc[d] = 1;
	}

	// preferred area for each face
	A0.resize(FNUM);

	// initialize nearest neighbor info
	NBX = -1;

	// seed random number generator
	srand48(seed);
}

// destructor
dpm::~dpm() {
	// clear all private vectors
	L.clear();
	pbc.clear();
	nv.clear();
	szList.clear();
	r.clear();
	x.clear();
	v.clear();
	F.clear();
	sb.clear();
	lb.clear();
	for (int i = 0; i < NBX; i++)
		nn.at(i).clear();
	nn.clear();
	head.clear();
	last.clear();
	list.clear();
    C_list.clear();

	if (posout.is_open())
		posout.close();
}

/******************************

	C E L L   S H A P E

	G E T T E R S

*******************************/

// get global vertex index gi given input cell index ci and vertex index vi
int dpm::gindex(int ci, int vi) {
	return szList[ci] + vi;
}

// get cell index ci and vertex index
void dpm::cindices(int &ci, int &vi, int gi) {
	for (int i = NCELLS - 1; i >= 0; i--) {
		if (gi >= szList[i]) {
			ci = i;
			vi = gi - szList[ci];
			break;
		}
	}
}

// get cell area
double dpm::volume(int ci) {
	// local variables
	int nf, ns1, ns2, ns3, gi, s_idx_start, ns_start;
	double V = 0.0;
    double             x1, y1, z1;
    double             x2, y2, z2;
    double             x3, y3, z3;
    double             dx12, dy12, dz12;
    double             dx13, dy13, dz13;
    double             Ax, Ay, Az;
    vector<double>     px(NVTOT, 0.0), py(NVTOT, 0.0), pz(NVTOT, 0.0);
    
    for (gi = 0; gi < NVTOT; gi++){
        ns_start = 3 * gi;
        px[gi] = x[ns_start];
        py[gi] = x[ns_start + 1];
        pz[gi] = x[ns_start + 2];
    }
    
    // loop over vertices of cell ci, get area by shoe-string method
    s_idx_start = szList.at(ci);
    for (nf = 0; nf < FNUM; nf++){
        ns1 = f_unit[nf][0] + s_idx_start;
        ns2 = f_unit[nf][1] + s_idx_start;
        ns3 = f_unit[nf][2] + s_idx_start;

        x1 = px[ns1]; y1 = py[ns1]; z1 = pz[ns1];
        x2 = px[ns2]; y2 = py[ns2]; z2 = pz[ns2];
        x3 = px[ns3]; y3 = py[ns3]; z3 = pz[ns3];

        dx12 = x2 - x1;
        dy12 = y2 - y1;
        dz12 = z2 - z1;
        dx13 = x3 - x1;
        dy13 = y3 - y1;
        dz13 = z3 - z1;

        Ax = (dy12 * dz13 - dz12 * dy13) * 0.5;
        Ay = (dz12 * dx13 - dx12 * dz13) * 0.5;
        Az = (dx12 * dy13 - dy12 * dx13) * 0.5;


        V += (Ax * x1 + Ay * y1 + Az * z1);
    }
    V /= 3.0;

	return abs(V);
}


// get cell perimeter
double dpm::area(int ci) {
	// local variables
    double A=0.0;
    double             x1, y1, z1;
    double             x2, y2, z2;
    double             x3, y3, z3;
    double             dx12, dy12, dz12;
    double             dx13, dy13, dz13;
    double             Ax, Ay, Az;
    int gi, ns_start, nf, ns1, ns2, ns3, s_idx_start;
    vector<double>     px(NVTOT, 0.0), py(NVTOT, 0.0), pz(NVTOT, 0.0);
    
    for (gi = 0; gi < NVTOT; gi++){
        ns_start = 3 * gi;
        px[gi] = x[ns_start];
        py[gi] = x[ns_start + 1];
        pz[gi] = x[ns_start + 2];
    }
    
    s_idx_start = szList.at(ci);
    for (nf = 0; nf < FNUM; nf++){
        ns1 = f_unit[nf][0] + s_idx_start;
        ns2 = f_unit[nf][1] + s_idx_start;
        ns3 = f_unit[nf][2] + s_idx_start;

        x1 = px[ns1]; y1 = py[ns1]; z1 = pz[ns1];
        x2 = px[ns2]; y2 = py[ns2]; z2 = pz[ns2];
        x3 = px[ns3]; y3 = py[ns3]; z3 = pz[ns3];

        dx12 = x2 - x1;
        dy12 = y2 - y1;
        dz12 = z2 - z1;
        dx13 = x3 - x1;
        dy13 = y3 - y1;
        dz13 = z3 - z1;

        Ax = (dy12 * dz13 - dz12 * dy13) * 0.5;
        Ay = (dz12 * dx13 - dx12 * dz13) * 0.5;
        Az = (dx12 * dy13 - dy12 * dx13) * 0.5;
        A += sqrt(Ax * Ax + Ay * Ay + Az * Az);
    }
	// return perimeter
	return A;
}


// get cell center of mass position
void dpm::com3D(int ci, double &cx, double &cy, double &cz) {
	// local variables
	int vi, gi, nvtmp;
    double xi, yi, zi, dx, dy, dz;

	// initial position: vi = 0
	nvtmp = nv.at(ci);
    gi = szList.at(ci);

    xi = x[NDIM * (gi)];
    yi = x[NDIM * (gi) + 1];
    zi = x[NDIM * (gi) + 2];
    
    cx = xi;
    cy = yi;
    cz = zi;
    
    for (vi = 1; vi < nvtmp; vi++){
        dx = x[NDIM * (gi + vi)] - xi;
        if (pbc[0])
            dx -= L[0] * round(dx / L[0]);
        xi += dx;
        cx += xi;
        
        dy = x[NDIM * (gi + vi) + 1] - yi;
        if (pbc[1])
            dy -= L[1] * round(dy / L[1]);
        yi += dy;
        cy += yi;
        
        dz = x[NDIM * (gi + vi) + 2] - zi;
        if (pbc[2])
            dz -= L[2] * round(dz / L[2]);
        zi += dz;
        cz += zi;
    }

	// take average to get com
	cx /= nvtmp;
	cy /= nvtmp;
    cz /= nvtmp;
    
    cx = fmod(cx,L[0]);
    cy = fmod(cy,L[1]);
    cz = fmod(cz,L[2]);
}

// get cell center of mass position
void dpm::cof3D(int ci, double &fx, double &fy, double &fz) {
    // local variables
    int vi, gi, nvtmp;
    
    // initial position: vi = 0
    nvtmp = nv.at(ci);
    gi = szList.at(ci);

    fx = F[NDIM * (gi)];
    fy = F[NDIM * (gi) + 1];
    fz = F[NDIM * (gi) + 2];
    
    for (vi = 1; vi < nvtmp; vi++){
        fx += F[NDIM * (gi + vi)];
        fy += F[NDIM * (gi + vi)+1];
        fz += F[NDIM * (gi + vi)+2];
    }
}

// get cell center of mass position
void dpm::cov3D(int ci, double &vx, double &vy, double &vz) {
    // local variables
    int vi, gi, nvtmp;

    // initial position: vi = 0
    nvtmp = nv.at(ci);
    gi = szList.at(ci);

    vx = v[NDIM * (gi)];
    vy = v[NDIM * (gi) + 1];
    vz = v[NDIM * (gi) + 2];
    
    for (vi = 1; vi < nvtmp; vi++){
        vx += v[NDIM * (gi + vi)];
        vy += v[NDIM * (gi + vi)+1];
        vz += v[NDIM * (gi + vi)+2];
    }
}

double dpm::NN3D(int gi){
    int gj, ci, vi, cj;
    double R, dx, dy, dz, x1, x2, y1, y2, z1, z2;
    double r_min=100, rj, r_max;
    r_max = 2.0 * (*max_element(r.begin(),r.end()));
    x1 = x[NDIM*gi];
    y1 = x[NDIM*gi+1];
    z1 = x[NDIM*gi+2];
    cindices(ci, vi, gi);
    for(gj=0;gj<NVTOT;gj++){
        cindices(cj, vi, gj);
        if (ci==cj) {
            continue;
        }
        x2 = x[NDIM*gj];
        y2 = x[NDIM*gj+1];
        z2 = x[NDIM*gj+2];
        
        dx = x1-x2;
        if (pbc[0])
            dx -= L[0] * round(dx / L[0]);
        if(dx < r_max && dx>0){
            dy = y1-y2;
            if (pbc[0])
                dy -= L[1] * round(dy / L[1]);
            if(dy<r_max && dy>0){
                dz = z1-z2;
                if (pbc[2])
                    dz -= L[2] * round(dz / L[2]);
                if(dz<r_max && dz>0){
                    R = sqrt(dx*dx+dy*dy+dz*dz)/(r[gi]+r[gj]) ;
                    if (r_min>R) {
                        r_min = R;
                        rj = r[gj];
                    }
                }
            }
        }
                
    }
    
    r_min = r_min/(r[gi]+rj) > x[gi*NDIM]/r[gi] ? x[gi*NDIM]/r[gi] : r_min/(r[gi]+rj);
    return r_min;
}

// get configuration packing fraction
double dpm::vertexPackingFraction3D() {
	int ci;
	double val, boxV, volumeSum = 0.0;

	// numerator
    for (ci = 0; ci < NCELLS; ci++){
        if (nv.at(ci)==1) {
            volumeSum += V0.at(ci);
        } else {
            volumeSum += volume(ci);
        }
    }

	// denominator
	boxV = L[0] * L[1] * L[2];

	// return packing fraction
	val = volumeSum / boxV;
	return val;
}

// get configuration "preferred" packing fraction
double dpm::vertexPreferredPackingFraction3D() {
    int ci;
    double val, boxV, volumeSum = 0.0;

    // numerator
    for (ci = 0; ci < NCELLS; ci++){
        if (nv.at(ci)==1) {
            volumeSum += V0.at(ci);
        } else {
            volumeSum += V0.at(ci);
        }
    }

    // denominator
    boxV = L[0] * L[1] * L[2];

    // return packing fraction
    val = volumeSum / boxV;
    return val;
}


/******************************

	I N I T I A L -

			I Z A T I O N

*******************************/

// initialize vertex shape parameters for SPECIFIC CELL ci
void dpm::initializeVertexShapeParameters(int ci, double calA0, double lenscale) {
	// local variables
	int gi, vi, nvtmp;

	// check that vertDOF has been assigned
	if (NVTOT <= 0){
		cerr << "	** ERROR: in initializeVertexShapeParameters, NVTOT not assigned. Ending here." << endl;
		exit(1);
	}
	if (vertDOF <= 0){
		cerr << "	** ERROR: in initializeVertexShapeParameters, vertDOF not assigned. Ending here." << endl;
		exit(1);
	}
	else if (nv.size() == 0){
		cerr << "	** ERROR: in initializeVertexShapeParameters, nv vector not assigned. Ending here." << endl;
    	exit(1);
	}

	// check that ci is within bounds
	if (ci >= NCELLS){
		cerr << "	** ERROR: in initializeVertexShapeParameters, ci = " << ci << ", but NCELLS = " << NCELLS << ". Will cause out of bounds error, so ending here. " << endl;
		exit(1);
	}

	// check that r have been resized correctly
	if (r.size() != NVTOT){
		cerr << "	** ERROR: in initializeVertexShapeParameters, r vector assigned incorrectly. r.size() = " << r.size() << ", while NVTOT = " << NVTOT << ". Ending here." << endl;
    	exit(1);
	}

	// since everything has been pre-allocated at this point, determine shape parameters based on input calA0
	nvtmp = nv.at(ci);

	V0.at(ci) = 1.0;

	// vertex radii
	gi = szList.at(ci);
    for (vi = 0; vi < nvtmp; vi++) {
        if (nvtmp == 1) {
            r.at(gi + vi) = pow(3.0/4.0/PI*V0.at(ci), 1.0/3.0) * lenscale;
            V0.at(ci) = 4.0/3.0*PI*r.at(gi + vi)*r.at(gi + vi)*r.at(gi + vi);
        } else {
            r.at(gi + vi) = D0_inter_unit * lenscale;
        }
    }
    
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
     */
}

// change size of particles
void dpm::scaleParticleSizes3D(double scaleFactor) {
	// local variables
	int ci, vi;
    double cx, cy, cz;
    double dx, dy, dz;
    double nvtmp;
    double s_idx;
    double scale_min = (L[0]-r.back())/L[0];
    scaleFactor = scaleFactor < scale_min ? scale_min : scaleFactor;
    L[0] *= scaleFactor;
    L[1] *= scaleFactor;
    L[2] *= scaleFactor;
    
    for(ci=0;ci<NCELLS;ci++){
        nvtmp = nv.at(ci);
        if (nvtmp==1) {
            x[ci*NDIM] *= scaleFactor;
            x[ci*NDIM + 1] *= scaleFactor;
            x[ci*NDIM + 2] *= scaleFactor;
        }
        else {
            com3D(ci, cx , cy, cz);
            dx = (scaleFactor -1) * cx;
            dy = (scaleFactor -1) * cy;
            dz = (scaleFactor -1) * cz;
            s_idx = szList[ci] * NDIM;
            for(vi=0;vi<nvtmp;vi++){
                x[s_idx + vi*NDIM] += dx;
                x[s_idx + vi*NDIM + 1] += dy;
                x[s_idx + vi*NDIM + 2] += dz;
            }
        }
    }
}

/******************************

	D P M  F O R C E 

			U P D A T E S

*******************************/

/******************************

	D P M  

		I N T E G R A T O R S

*******************************/

void dpm::setdt(double dt0) {
	dt = dt0;
}
