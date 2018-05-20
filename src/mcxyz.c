/********************************************
 *
 *  mcxyz.c,	in the C programming language
 *
 * created 2010, 2012 by
 *	Steven L. JACQUES
 *  Ting LI
 *	Oregon Health & Science University
 *
 * Log:
 *  2010: Written by Ting Li based on Steve's mcsub.c.
 *  2010-12-30: Use Steve Jacques' FindVoxelFace().
 *  2012-05-08: Reorganized by Steve Jacques.
 *  2014-01-09: Edited to included a uniformly distributed light source over the entire surface by Mathias Christensen, DTU Fotonik.
 *  2017-04-20: Overhauled by Anders K. Hansen and Dominik Marti, DTU Fotonik. Fundamental method remained unchanged.
 *      Uses the Mersenne Twister for random number generation.
 *  2017-06-07: Adapted to MATLAB mex file generation by Anders K. Hansen, DTU Fotonik
 *  2018-04-19: Added support for illuminating with an isotropically emitting 3D distribution of sources, such as a collection of fluorescing emitters
 *
 ** COMPILING ON WINDOWS
 * Can be compiled in MATLAB with "mex COPTIMFLAGS='$COPTIMFLAGS -Ofast -fopenmp -std=c11 -Wall -pedantic' LDOPTIMFLAGS='$LDOPTIMFLAGS -Ofast -fopenmp -std=c11 -Wall -pedantic' -outdir private .\src\mcxyz.c ".\src\libut.lib""
 *
 * To get the MATLAB C compiler to work, try this:
 * 1. Go to MATLAB's addon manager and tell it to install the "Support for MinGW-w64 compiler"
 * 2. Type "mex -setup" in the MATLAB command window and ensure that MATLAB has set the C compiler to MinGW64
 * 3. If you're using an older MATLAB version, you may need to copy the files "libgomp.a" and "libgomp.spec" to the folder with a path similar to "C:\ProgramData\MATLAB\SupportPackages\R2017a\MW_MinGW_4_9\lib\gcc\x86_64-w64-mingw32\4.9.2"
 * 4. mex should now be able to compile the code using the above command but in order to run, it needs to have the file "libgomp_64-1.dll" copied to the same folder as the mex file.
 *
 ** COMPILING ON MAC
 * As of June 2017, the macOS compiler doesn't support libut (for ctrl+c 
 * breaking) or openmp (for multithreading).
 * Compile in MATLAB with "mex COPTIMFLAGS='$COPTIMFLAGS -Ofast -std=c11 -Wall -pedantic' LDOPTIMFLAGS='$LDOPTIMFLAGS -Ofast -std=c11 -Wall -pedantic' -outdir private ./src/mcxyz.c"
 *
 * To get the MATLAB C compiler to work, try this:
 * 1. Install XCode from the App Store
 * 2. Type "mex -setup" in the MATLAB command window
 ********************************************/

#include "mex.h"
#include <math.h>
#include <time.h>
#ifdef _WIN32 // This is defined on both win32 and win64 systems. We use this preprocessor condition to avoid loading openmp or libut on, e.g., Mac
#include <omp.h>
#endif

#define DSFMT_MEXP 19937 // Mersenne exponent for dSFMT
#include "dSFMT-src-2.2.3/dSFMT.c" // Double precision SIMD oriented Fast Mersenne Twister(dSFMT)

#include "lambert.c" // For calculating the Lambert W function, originally part of the GNU Scientific Library, created by K. Briggs, G. Jungman and B. Gough and slightly modified by A. Hansen for easier mcxyz integration

#define PI          3.14159265358979323846
#define THRESHOLD   0.01		// used in roulette
#define CHANCE      0.1  		// used in roulette
#define SIGN(x)     ((x)>=0 ? 1:-1)
#define RandomNum   dsfmt_genrand_open_close(&dsfmt) // Calls for a random number in (0,1]

// Function for rotating a point at coordinates (x,y,z) an angle theta around an axis with direction unit vector (ux,uy,uz).
void axisrotate(double *x, double *y, double *z, double ux, double uy, double uz, double theta) {	
    double xin = *x, yin = *y, zin = *z;
	double st = sin(theta), ct = cos(theta);
	
	*x = (ux*ux*(1-ct) +    ct)*xin	 +  (ux*uy*(1-ct) - uz*st)*yin  +  (ux*uz*(1-ct) + uy*st)*zin;
	*y = (uy*ux*(1-ct) + uz*st)*xin	 +  (uy*uy*(1-ct) +    ct)*yin  +  (uy*uz*(1-ct) - ux*st)*zin;
	*z = (uz*ux*(1-ct) - uy*st)*xin	 +  (uz*uy*(1-ct) + ux*st)*yin  +  (uz*uz*(1-ct) +    ct)*zin;
}

#ifdef _WIN32
extern bool utIsInterruptPending(); // Allows catching ctrl+c while executing the mex function
#endif

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, mxArray const *prhs[]) {
    // Extracting quantities from MCinput struct
    double  simulationTimeRequested     = *mxGetPr(mxGetField(prhs[0],0,"simulationTime"));
    int     beamtypeFlag                = *mxGetPr(mxGetField(prhs[0],0,"beamtypeFlag"));
    int     boundaryFlag                = *mxGetPr(mxGetField(prhs[0],0,"boundaryFlag"));
    double  xFocus                      = *mxGetPr(mxGetField(prhs[0],0,"xFocus"));
    double  yFocus                      = *mxGetPr(mxGetField(prhs[0],0,"yFocus"));
    double  zFocus                      = *mxGetPr(mxGetField(prhs[0],0,"zFocus"));
    double  ux0                         = *mxGetPr(mxGetField(prhs[0],0,"ux0"));
    double  uy0                         = *mxGetPr(mxGetField(prhs[0],0,"uy0"));
    double  uz0                         = *mxGetPr(mxGetField(prhs[0],0,"uz0"));
    double  waist                       = *mxGetPr(mxGetField(prhs[0],0,"waist"));
    double  divergence                  = *mxGetPr(mxGetField(prhs[0],0,"divergence"));
    double  dx                          = *mxGetPr(mxGetField(prhs[0],0,"dx"));
    double  dy                          = *mxGetPr(mxGetField(prhs[0],0,"dy"));
    double  dz                          = *mxGetPr(mxGetField(prhs[0],0,"dz"));
    
    mxArray *sourceDist = mxGetField(prhs[0],0,"sourceDistribution");
    double *S           = mxGetData(sourceDist);
    double sourcesum    = 0;
    double *S_CDF       = 0;
    
    mxArray *tissueList = mxGetField(prhs[0],0,"tissueList");
    mxArray *T          = mxGetField(prhs[0],0,"T");
    mwSize const *dimPtr = mxGetDimensions(T);
    int     nx = dimPtr[0];
    int     ny = dimPtr[1];
    int     nz = dimPtr[2];
    int     nT = mxGetN(tissueList);
    unsigned char *v = mxGetData(T);
    
    long i;                      // General-purpose non-thread-specific index variable
    double 	muav[nT];            // muav[0:nT-1], absorption coefficient of ith tissue type
	double 	musv[nT];            // scattering coeff. 
	double 	gv[nT];              // anisotropy of scattering
    for(i=0;i<nT;i++) {
        muav[i] = *mxGetPr(mxGetField(tissueList,i,"mua"));
        musv[i] = *mxGetPr(mxGetField(tissueList,i,"mus"));
        gv[i]   = *mxGetPr(mxGetField(tissueList,i,"g"));
    }
    
    // Prepare output matrix
    plhs[0] = mxCreateNumericArray(3,dimPtr,mxDOUBLE_CLASS,mxREAL);
    double  *F = mxGetPr(plhs[0]);
    for(i=0;i<nx*ny*nz;i++) F[i] = 0;
    
	// other variables
	long long	nPhotons;                   // number of photons in simulation
	int     pctProgress = 0;                // Simulation progress in percent
    int     newPctProgress;    
    bool    ctrlc_caught = false;           // Has a ctrl+c been passed from MATLAB?
	double  photonKillRange = 6;
    /* The photonKillRange determines the size of the region that 
     * photons are launched in for the infinite plane wave without 
     * boundaries but also the region that photons are allowed to stay 
     * alive in (if outside, the probability of returning to the region
     * of interest is judged as too low) */
    
	// launch parameters
	double vx0 = 1, vy0 = 0, vz0 = 0;          // normal vector to photon trajectory

    // time
	double	simulationTimeSpent;               // Spent time durations of computation.
	struct  timespec simulationTimeStart, simulationTimeCurrent;
	
    // Display parameters
    printf("---------------------------------------------------------\n");
    printf("(nx,ny,nz) = (%d,%d,%d), (dx,dy,dz) = (%0.1e,%0.1e,%0.1e) cm\n",nx,ny,nz,dx,dy,dz);

    if(sourceDist) printf("Using isotropically emitting 3D distribution as light source\n");
    else {
        switch(beamtypeFlag) {
            case 0:
                printf("Using top-hat focus, top-hat far field beam\n");
                break;
            case 1:
                printf("Using Gaussian focus, Gaussian far field beam\n");
                break;
            case 2:
                printf("Using isotropic point source\n");
                break;
            case 3:
                printf("Using infinite plane wave\n");
                break;
            case 4:
                printf("Using pencil beam\n");
                break;
            case 5:
                printf("Using top-hat focus, Gaussian far field beam\n");
                break;
            case 6:
                printf("Using Gaussian focus, top-hat far field beam\n");
                break;
            case 7:
                printf("Using Laguerre-Gaussian LG01 focus, LG01 far field beam\n");
                break;
        }

        printf("Focus location = (%0.1e,%0.1e,%0.1e) cm\n",xFocus,yFocus,zFocus);
        if(beamtypeFlag!=2) printf("Beam direction = (%0.3f,%0.3f,%0.3f)\n",ux0,uy0,uz0);
        if(beamtypeFlag!=2 && beamtypeFlag!=3 && beamtypeFlag!=4) {
            printf("Beam waist = %0.1e cm, divergence = %0.1e rad\n",waist,divergence);
        }
    }
    
    if(boundaryFlag==0) printf("boundaryFlag = 0, so no boundaries\n");
    else if(boundaryFlag==1) printf("boundaryFlag = 1, so escape at all boundaries\n");    
	else if(boundaryFlag==2) printf("boundaryFlag = 2, so escape at surface only\n");    
	else {
        printf("Error: Invalid boundaryFlag\n");
        return;
    }
    
	// INITIALIZATIONS
	if(sourceDist) {
        S_CDF = mxMalloc((nx*ny*nz + 1)*sizeof(double));
        S_CDF[0] = 0;
        for(i=1;i<(nx*ny*nz+1);i++) S_CDF[i] = S_CDF[i-1] + S[i-1];
        sourcesum = S_CDF[nx*ny*nz];
        for(i=1;i<(nx*ny*nz+1);i++) S_CDF[i] /= sourcesum;
    } else if(uz0!=1) {
		vx0 = -uy0/sqrt(uy0*uy0 + ux0*ux0);
		vy0 = ux0/sqrt(uy0*uy0 + ux0*ux0);
		vz0 = 0;
	}

    printf("Simulation duration = %0.2f min\n\n",simulationTimeRequested);
	
	// RUN - Launch photons, initializing each one before propagation
	printf("Calculating...   0%% done");
    mexEvalString("drawnow;");
	
    nPhotons = 0;
	
    // ============================ MAJOR CYCLE ========================
	clock_gettime(CLOCK_MONOTONIC, &simulationTimeStart);
	
    #ifdef _WIN32
    #pragma omp parallel num_threads(omp_get_num_procs())
    #endif
	{
		double	ux, uy, uz;         // photon trajectory as unit vector composants
        double  rand;               // Random number used for lookup into CDF of 3D source
        long    d;                  // Used in CDF lookup
		double  ux_temp, uy_temp;   // temporary values used during SPIN
		double	s;                  // step sizes. s = -log(RND)/mus [cm]
		double  stepLeft;           // dimensionless step size
		double	costheta;           // cos(theta)
		double  sintheta;           // sin(theta)
		double	cospsi;             // cos(psi)
		double  sinpsi;             // sin(psi)
		double	psi;                // azimuthal angle
		double	photonWeight;       // photon weight
		double	absorb;             // weighted deposited in a step due to absorption
		bool    photonAlive;        // flag, true or false
		bool    sameVoxel;          // Are they in the same voxel?
		double	mua;                // absorption coefficient [cm^-1]
		double	mus;                // scattering coefficient [cm^-1]
		double	g;                  // anisotropy [-]
		double  xTarget, yTarget, zTarget;
		double	wx0, wy0, wz0;      // normal vector 2 to photon trajectory. Not necessarily normal to v0.
		double	r, theta, phi;
        long    j;                  // General-purpose thread-specific index variable
        double  sx, sy, sz, sxyzmin;// Distances along current trajectory to next voxel boundary plane parallel to the yz, xz and xy planes, respectively, and the minimum of the three.
		double 	ix=0, iy=0, iz=0;   // Fractional x, y, or z indices. Used to track photons
		bool    photonInsideVolume; // Is the photon inside the cuboid?
		dsfmt_t dsfmt;			    // Thread-specific "state" of dSFMT pseudo-random number generator

        #ifdef _WIN32
		dsfmt_init_gen_rand(&dsfmt,simulationTimeStart.tv_nsec + omp_get_thread_num()); // Seed random number generator
        #else
		dsfmt_init_gen_rand(&dsfmt,simulationTimeStart.tv_nsec); // Seed random number generator
        #endif
		while((pctProgress < 100) && !ctrlc_caught) {
			// LAUNCH - Initialize photon position and trajectory

			photonWeight = 1;
			photonAlive = true;
			sameVoxel = false; // Photon is initialized as if it has just entered the voxel it's created in, to ensure proper initialization of voxel index and ray properties
			
            #ifdef _WIN32
			#pragma omp atomic
            #endif
            nPhotons += 1;
			
            if(sourceDist) { // If a 3D source distribution was passed in as part of the input struct (used for fluorescence)...
                // ... then search the cumulative distribution function via binary tree method to find the voxel to start the photon in
                rand = RandomNum;
                d = nx*ny*nz-1; // Number of elements in the current section to be searched, minus one
                j = d/2; // Index of the middle element of the current section
                while(!(S_CDF[j] < rand && S_CDF[j+1] >= rand)) { // Binary tree search
                    if(S_CDF[j] >= rand) {
                        j = j + (d-2)/4 - d/2;
                        d = d/2 - 1;
                    } else {
                        d = d - d/2 - 1;
                        j = j + d/2 + 1;
                    }
                }
                
                ix = j%nx    + 1 - RandomNum;
                iy = j/nx%ny + 1 - RandomNum;
                iz = j/nx/ny + 1 - RandomNum;
                
                costheta = 1 - 2*RandomNum;
                sintheta = sqrt(1 - costheta*costheta);
                psi = 2*PI*RandomNum;
                
                ux = sintheta*cos(psi);
                uy = sintheta*sin(psi);
                uz = costheta;
            } else switch (beamtypeFlag) {
                case 0: // top-hat focus, top-hat far field beam
                    r		= waist*sqrt(RandomNum); // for target calculation
                    phi		= RandomNum*2*PI;
                    wx0		= vx0;
                    wy0		= vy0;
                    wz0		= vz0;
                    axisrotate(&wx0,&wy0,&wz0,ux0,uy0,uz0,phi); // w0 unit vector now points in the direction from focus center point to ray target point

                    xTarget	= xFocus + r*wx0;
                    yTarget	= yFocus + r*wy0;
                    zTarget = zFocus + r*wz0;

                    theta	= divergence*sqrt(RandomNum); // for trajectory calculation. The sqrt is valid within paraxial approximation.
                    phi		= RandomNum*2*PI;
                    wx0		= vx0;
                    wy0		= vy0;
                    wz0		= vz0;
                    axisrotate(&wx0,&wy0,&wz0,ux0,uy0,uz0,phi); // w0 unit vector is now normal to both beam center axis and to ray propagation direction. Angle from v0 to w0 is phi.

                    ux		= ux0;
                    uy		= uy0;
                    uz		= uz0;
                    axisrotate(&ux,&uy,&uz,wx0,wy0,wz0,theta); // ray propagation direction is found by rotating beam center axis an angle theta around w0

                    ix		= (xTarget - zTarget*ux/uz)/dx + nx/2.0; // the coordinates for the ray starting point is the intersection of the ray with the z = 0 surface
                    iy		= (yTarget - zTarget*uy/uz)/dy + ny/2.0;
                    iz		= 0;
                break;
                case 1: // Gaussian focus, Gaussian far field beam
                    r		= waist*sqrt(-0.5*log(RandomNum)); // for target calculation
                    phi		= RandomNum*2*PI;
                    wx0		= vx0;
                    wy0		= vy0;
                    wz0		= vz0;
                    axisrotate(&wx0,&wy0,&wz0,ux0,uy0,uz0,phi); // w0 unit vector now points in the direction from focus center point to ray target point

                    xTarget	= xFocus + r*wx0;
                    yTarget	= yFocus + r*wy0;
                    zTarget = zFocus + r*wz0;

                    theta	= divergence*sqrt(-0.5*log(RandomNum)); // for trajectory calculation. The sqrt is valid within paraxial approximation.
                    phi		= RandomNum*2*PI;
                    wx0		= vx0;
                    wy0		= vy0;
                    wz0		= vz0;
                    axisrotate(&wx0,&wy0,&wz0,ux0,uy0,uz0,phi); // w0 unit vector is now normal to both beam center axis and to ray propagation direction. Angle from v0 to w0 is phi.

                    ux		= ux0;
                    uy		= uy0;
                    uz		= uz0;
                    axisrotate(&ux,&uy,&uz,wx0,wy0,wz0,theta); // ray propagation direction is found by rotating beam center axis an angle theta around w0

                    ix		= (xTarget - zTarget*ux/uz)/dx + nx/2.0; // the coordinates for the ray starting point is the intersection of the ray with the z = 0 surface
                    iy		= (yTarget - zTarget*uy/uz)/dy + ny/2.0;
                    iz		= 0;
                break;
                case 2: // isotropically emitting point source
                    ix = xFocus/dx + nx/2.0;
                    iy = yFocus/dy + ny/2.0;
                    iz = zFocus/dz;
                    
                    costheta = 1 - 2*RandomNum;
                    sintheta = sqrt(1 - costheta*costheta);
                    psi = 2*PI*RandomNum;
                    
                    ux = sintheta*cos(psi);
                    uy = sintheta*sin(psi);
                    uz = costheta;
                break;
                case 3: // infinite plane wave
                    if(boundaryFlag==1) {
                        ix = nx*(RandomNum-0.5) + nx/2.0; // Generates a random ix coordinate within the cuboid
                        iy = ny*(RandomNum-0.5) + ny/2.0; // Generates a random iy coordinate within the cuboid
                    }
                    else {
                        ix = photonKillRange*nx*(RandomNum-0.5) + nx/2.0; // Generates a random ix coordinate within an interval photonKillRange times the cuboid's size
                        iy = photonKillRange*ny*(RandomNum-0.5) + ny/2.0; // Generates a random iy coordinate within an interval photonKillRange times the cuboid's size
                    }
                    iz = 0;
                    ux = ux0;
                    uy = uy0;
                    uz = uz0;
                break;
                case 4: // pencil beam
                    ix		= (xFocus - zFocus*ux/uz)/dx + nx/2.0;
                    iy		= (yFocus - zFocus*uy/uz)/dy + ny/2.0;
                    iz		= 0;
                    ux	= ux0;
                    uy	= uy0;
                    uz	= uz0;
                break;
                case 5: // top-hat focus, Gaussian far field beam
                    r		= waist*sqrt(RandomNum); // for target calculation
                    phi		= RandomNum*2*PI;
                    wx0		= vx0;
                    wy0		= vy0;
                    wz0		= vz0;
                    axisrotate(&wx0,&wy0,&wz0,ux0,uy0,uz0,phi); // w0 unit vector now points in the direction from focus center point to ray target point

                    xTarget	= xFocus + r*wx0;
                    yTarget	= yFocus + r*wy0;
                    zTarget = zFocus + r*wz0;

                    theta	= divergence*sqrt(-0.5*log(RandomNum)); // for trajectory calculation. The sqrt is valid within paraxial approximation.
                    phi		= RandomNum*2*PI;
                    wx0		= vx0;
                    wy0		= vy0;
                    wz0		= vz0;
                    axisrotate(&wx0,&wy0,&wz0,ux0,uy0,uz0,phi); // w0 unit vector is now normal to both beam center axis and to ray propagation direction. Angle from v0 to w0 is phi.

                    ux		= ux0;
                    uy		= uy0;
                    uz		= uz0;
                    axisrotate(&ux,&uy,&uz,wx0,wy0,wz0,theta); // ray propagation direction is found by rotating beam center axis an angle theta around w0

                    ix		= (xTarget - zTarget*ux/uz)/dx + nx/2.0; // the coordinates for the ray starting point is the intersection of the ray with the z = 0 surface
                    iy		= (yTarget - zTarget*uy/uz)/dy + ny/2.0;
                    iz		= 0;
                break;
                case 6: // Gaussian focus, top-hat far field beam
                    r		= waist*sqrt(-0.5*log(RandomNum)); // for target calculation
                    phi		= RandomNum*2*PI;
                    wx0		= vx0;
                    wy0		= vy0;
                    wz0		= vz0;
                    axisrotate(&wx0,&wy0,&wz0,ux0,uy0,uz0,phi); // w0 unit vector now points in the direction from focus center point to ray target point

                    xTarget	= xFocus + r*wx0;
                    yTarget	= yFocus + r*wy0;
                    zTarget = zFocus + r*wz0;

                    theta	= divergence*sqrt(RandomNum); // for trajectory calculation. The sqrt is valid within paraxial approximation.
                    phi		= RandomNum*2*PI;
                    wx0		= vx0;
                    wy0		= vy0;
                    wz0		= vz0;
                    axisrotate(&wx0,&wy0,&wz0,ux0,uy0,uz0,phi); // w0 unit vector is now normal to both beam center axis and to ray propagation direction. Angle from v0 to w0 is phi.

                    ux		= ux0;
                    uy		= uy0;
                    uz		= uz0;
                    axisrotate(&ux,&uy,&uz,wx0,wy0,wz0,theta); // ray propagation direction is found by rotating beam center axis an angle theta around w0

                    ix		= (xTarget - zTarget*ux/uz)/dx + nx/2.0; // the coordinates for the ray starting point is the intersection of the ray with the z = 0 surface
                    iy		= (yTarget - zTarget*uy/uz)/dy + ny/2.0;
                    iz		= 0;
                break;
                case 7: // Laguerre-Gaussian LG01 focus, LG01 far field beam
                    r		= waist*sqrt((gsl_sf_lambert_Wm1(-RandomNum*exp(-1))+1)/(-2))/1.50087; // for target calculation
                    phi		= RandomNum*2*PI;
                    wx0		= vx0;
                    wy0		= vy0;
                    wz0		= vz0;
                    axisrotate(&wx0,&wy0,&wz0,ux0,uy0,uz0,phi); // w0 unit vector now points in the direction from focus center point to ray target point

                    xTarget	= xFocus + r*wx0;
                    yTarget	= yFocus + r*wy0;
                    zTarget = zFocus + r*wz0;

                    theta	= divergence*sqrt((gsl_sf_lambert_Wm1(-RandomNum*exp(-1))+1)/(-2))/1.50087; // for trajectory calculation
                    ux		= ux0;
                    uy		= uy0;
                    uz		= uz0;
                    axisrotate(&ux,&uy,&uz,wx0,wy0,wz0,theta); // ray propagation direction is found by rotating beam center axis an angle theta around w0

                    ix		= (xTarget - zTarget*ux/uz)/dx + nx/2.0; // the coordinates for the ray starting point is the intersection of the ray with the z = 0 surface
                    iy		= (yTarget - zTarget*uy/uz)/dy + ny/2.0;
                    iz		= 0;
                break;
            }
            
			// HOP_DROP_SPIN_CHECK - Propagate one photon until it dies as determined by ROULETTE
			while(photonAlive) {
                // Calculate distances to next voxel boundary planes
                sx = ux? (floor(ix) + (ux>0) - ix)*dx/ux : INFINITY;
                sy = uy? (floor(iy) + (uy>0) - iy)*dy/uy : INFINITY;
                sz = uz? (floor(iz) + (uz>0) - iz)*dz/uz : INFINITY;
                
				// HOP - Take step to new position
				stepLeft	= -log(RandomNum); // dimensionless step
				
				do { // while(stepLeft>0)
					if(!sameVoxel) { // If photon has just entered a new voxel or has just been launched
                        photonInsideVolume = (ix < nx && ix >= 0 && iy < ny && iy >= 0 && iz < nz && iz >= 0);
                        
                        switch (boundaryFlag) {
                            case 0:
                                photonAlive = (fabs(ix/nx - 1.0/2) <  photonKillRange/2 &&
                                               fabs(iy/ny - 1.0/2) <  photonKillRange/2 &&
                                               fabs(iz/nz - 1.0/2) <  photonKillRange/2);
                                break;
                            case 1:
                                photonAlive = photonInsideVolume;
                                break;
                            case 2:
                                photonAlive = (fabs(ix/nx - 1.0/2) <  photonKillRange/2 &&
                                               fabs(iy/ny - 1.0/2) <  photonKillRange/2 &&
                                                    iz/nz - 1.0/2  <  photonKillRange/2 &&
                                                    iz             >= 0);
                                break;
                        }
						if(!photonAlive) break;
						
						/* Get tissue voxel properties of current position.
						 * If photon is outside cuboid, properties are those of
                         * the closest defined voxel. */
                        j = ((iz < 0)? 0: ((iz >= nz)? nz-1: floor(iz)))*nx*ny +
                            ((iy < 0)? 0: ((iy >= ny)? ny-1: floor(iy)))*nx    +
                            ((ix < 0)? 0: ((ix >= nx)? nx-1: floor(ix))); // Index values are restrained to integers in the interval [0,n-1]
						mua 	= muav[v[j]];
						mus 	= musv[v[j]];
						g 		= gv[v[j]];
					}
					
					s     = stepLeft/mus; // Step size [cm]
                    
                    sxyzmin = fmin(sx,fmin(sy,sz));
                    if(s < sxyzmin) { // s is smallest, so the photon stays within the voxel
                        sameVoxel = true;
                        stepLeft = 0;
                    } else if(s == sxyzmin) { // photon gets scattered exactly on a voxel boundary
                        sameVoxel = false;
                        stepLeft = 0;
                    } else { // photon crosses voxel boundary
                        sameVoxel = false;
                        s = sxyzmin;
                        stepLeft -= s*mus;
                    }

                    ix = (sx==s)? ((ux>0)? floor(ix+1): nextafter(floor(ix),-INFINITY)): (ix + s*ux/dx); // The nextafter function is used to move the photon the smallest possible amount (machine precision) into the voxel.
                    iy = (sy==s)? ((uy>0)? floor(iy+1): nextafter(floor(iy),-INFINITY)): (iy + s*uy/dy);
                    iz = (sz==s)? ((uz>0)? floor(iz+1): nextafter(floor(iz),-INFINITY)): (iz + s*uz/dz);
                    
                    sx = (sx==s)? dx/fabs(ux): sx - s;
                    sy = (sy==s)? dy/fabs(uy): sy - s;
                    sz = (sz==s)? dz/fabs(uz): sz - s;
                    
                    // DROP - Drop photon weight into local bin
                    absorb = photonWeight*(1 - exp(-mua*s));   // photon weight absorbed at this step
                    photonWeight -= absorb;					   // decrement WEIGHT by amount absorbed
                    if(photonInsideVolume) {	// only save data if the photon is inside simulation cuboid
                        #ifdef _WIN32
                        #pragma omp atomic
                        #endif
                            F[j] += absorb;
                    }
				} while(stepLeft>0);
				
				/**** CHECK ROULETTE 
				 If photon weight below THRESHOLD, then terminate photon using Roulette technique.
				 Photon has CHANCE probability of having its weight increased by factor of 1/CHANCE,
				 and 1-CHANCE probability of terminating. */
				if(photonWeight < THRESHOLD) {
					if(RandomNum <= CHANCE) photonWeight /= CHANCE;
					else photonAlive = false;
				}
				
				if(photonAlive) {
					// SPIN - Scatter photon into new trajectory
					// Sample for costheta using Henyey-Greenstein scattering
                    costheta = g? (1 + g*g - pow((1 - g*g)/(1 - g + 2*g*RandomNum),2))/(2*g) : 2*RandomNum - 1;
					sintheta = sqrt(1 - costheta*costheta);
					
					psi = 2*PI*RandomNum;
					cospsi = cos(psi);
					sinpsi = sin(psi);
                    
                    if(fabs(uz) < 1) {
						ux_temp =  sintheta*(ux*uz*cospsi - uy*sinpsi)/sqrt(ux*ux + uy*uy) + ux*costheta;
						uy_temp =  sintheta*(uy*uz*cospsi + ux*sinpsi)/sqrt(ux*ux + uy*uy) + uy*costheta;
						uz      = -sintheta*(      cospsi            )*sqrt(ux*ux + uy*uy) + uz*costheta;
                        uy      = uy_temp;
                        ux      = ux_temp;
                    } else {
						ux = sintheta*cospsi;
						uy = sintheta*sinpsi;
						uz = costheta*SIGN(uz);
                    }
				}
			} // while(photonAlive)
			
            #ifdef _WIN32
			#pragma omp master
            #endif
			{
				// Check whether ctrl+c has been pressed
                #ifdef _WIN32
                if(utIsInterruptPending()) {
                    ctrlc_caught = true;
                    printf("\nCtrl+C detected, stopping.");
                }
                #endif
                clock_gettime(CLOCK_MONOTONIC, &simulationTimeCurrent);
 
				// Print out message about progress.
                newPctProgress = 100*(simulationTimeCurrent.tv_sec  - simulationTimeStart.tv_sec +
                                     (simulationTimeCurrent.tv_nsec - simulationTimeStart.tv_nsec)/1e9) / (simulationTimeRequested*60);
				if((newPctProgress != pctProgress) && !ctrlc_caught) {
					pctProgress = newPctProgress;
                    printf("\b\b\b\b\b\b\b\b\b%3.i%% done", pctProgress);
                    mexEvalString("drawnow;");
				}
			}
		} // while((pctProgress < 100) && !ctrlc_caught)
	}
    
	clock_gettime(CLOCK_MONOTONIC, &simulationTimeCurrent);
	simulationTimeSpent = simulationTimeCurrent.tv_sec  - simulationTimeStart.tv_sec +
                         (simulationTimeCurrent.tv_nsec - simulationTimeStart.tv_nsec)/1e9;
	printf("\nSimulated %0.1e photons at a rate of %0.1e photons per minute\n",(double)nPhotons, nPhotons*60/simulationTimeSpent);
    printf("---------------------------------------------------------\n");
    mexEvalString("drawnow;");
    
    if(sourceDist) mxFree(S_CDF);
    
    // SAVE - Convert data to relative fluence rate [cm^-2] and return
    
    // Normalize deposition to yield either absolute or relative fluence rate (F).
    if(sourceDist) { // For a 3D source distribution (e.g., fluorescence)
		for(i=0; i<nx*ny*nz;i++) F[i] /= nPhotons*muav[v[i]]/sourcesum; // Absolute fluence rate [W/cm^2]
    } else if(beamtypeFlag == 3 && boundaryFlag != 1) { // For infinite plane wave launched into volume without absorbing walls
		for(i=0; i<nx*ny*nz;i++) F[i] /= dx*dy*dz*nPhotons/(photonKillRange*photonKillRange)*muav[v[i]]; // Normalized fluence rate [W/cm^2/W.incident]
	} else {
		for(i=0; i<nx*ny*nz;i++) F[i] /= dx*dy*dz*nPhotons*muav[v[i]]; // Normalized fluence rate [W/cm^2/W.incident]
	}
	
    return;
}
