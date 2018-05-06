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
 * 3. Copy the files "libgomp.a" and "libgomp.spec" to the folder with a path similar to "C:\ProgramData\MATLAB\SupportPackages\R2017a\MW_MinGW_4_9\lib\gcc\x86_64-w64-mingw32\4.9.2"
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

#include "lambert.c" // For calculating the Lambert W function, originally part of the GNU Scientific Library, modified by Anders for easier mcxyz integration

#define ls          1.0E-7      /* Moving photon a little bit off the voxel face */
#define	PI          3.1415926
#define THRESHOLD   0.01		/* used in roulette */
#define CHANCE      0.1  		/* used in roulette */
#define SIGN(x)     ((x)>=0 ? 1:-1)
#define RandomNum   dsfmt_genrand_open_close(&dsfmt) /* Calls for a random number in (0,1] */
#define ONE_MINUS_COSZERO 1.0E-12   /* If 1-cos(theta) <= ONE_MINUS_COSZERO, fabs(theta) <= 1e-6 rad. */
/* If 1+cos(theta) <= ONE_MINUS_COSZERO, fabs(PI-theta) <= 1e-6 rad. */

// Function for rotating a point at coordinates (x,y,z) an angle theta around an axis with direction unit vector (ux,uy,uz).
void axisrotate(double *x, double *y, double *z, double ux, double uy, double uz, double theta) {	
    double xin = *x, yin = *y, zin = *z;
	double st = sin(theta), ct = cos(theta);
	
	*x = (ct + ux*ux*(1-ct))   *xin	 +  (ux*uy*(1-ct) - uz*st)*yin  +  (ux*uz*(1-ct) + uy*st)*zin;
	*y = (uy*ux*(1-ct) + uz*st)*xin	 +  (ct + uy*uy*(1-ct))   *yin  +  (uy*uz*(1-ct) - ux*st)*zin;
	*z = (uz*ux*(1-ct) - uy*st)*xin	 +  (uz*uy*(1-ct) + ux*st)*yin  +  (ct + uz*uz*(1-ct))   *zin;
}

// Distance to next voxel face along current trajectory
double FindVoxelFace(double x1, double y1, double z1, double dx, double dy, double dz, int nx, int ny, double ux, double uy, double uz) {	
    int ix = floor(nx/2.0 + x1/dx) + (ux>=0);
    int iy = floor(ny/2.0 + y1/dy) + (uy>=0);
    int iz = floor(         z1/dz) + (uz>=0);

    double xs = ux? ((ix - nx/2.0)*dx - x1)/ux : INFINITY;
    double ys = uy? ((iy - ny/2.0)*dy - y1)/uy : INFINITY;
    double zs = uz? ((iz         )*dz - z1)/uz : INFINITY;
    
    return fmin(xs,fmin(ys,zs));
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
    
    long i;        /* General-purpose non-thread-specific index variable */
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
    for(i=0;i<nx*ny*nz;i++) {
        F[i] = 0;
    }
    
	/* other variables */
	long long	nPhotons;                   /* number of photons in simulation */
	int     pctProgress = 0;                /* Simulation progress in percent */
    int     newPctProgress;    
    bool    ctrlc_caught = false;           // Has a ctrl+c been passed from MATLAB?
	double  infWaveScaleFactor = 6;
    /* The infWaveScaleFactor determines the size of the region that 
     * photons are launched in for the infinite plane wave without 
     * boundaries but also the region that photons are allowed to stay 
     * alive in (if outside, the probability of returning to the region
     * of interest is judged as too low)
     */
    
	/* launch parameters */
	double vx0 = 1, vy0 = 0, vz0 = 0;  /* normal vector to photon trajectory */

    /* time */
	double	simulationTimeSpent;               // Spent time durations of computation.
	time_t	currentCalendarTime;
	struct  timespec simulationTimeStart, simulationTimeCurrent;
	
    // Display parameters
    printf("nx = %d, dx = %0.8f [cm]\n",nx,dx);
    printf("ny = %d, dy = %0.8f [cm]\n",ny,dy);
    printf("nz = %d, dz = %0.8f [cm]\n\n",nz,dz);

    if (sourceDist) printf("Using isotropically emitting 3D distribution as light source\n\n");
    else {
        printf("beamtypeFlag = %d\n",beamtypeFlag);
        switch(beamtypeFlag) {
            case 0:
                printf("Using top-hat focus, top-hat far field beam\n\n");
                break;
            case 1:
                printf("Using Gaussian focus, Gaussian far field beam\n\n");
                break;
            case 2:
                printf("Using isotropic point source\n\n");
                break;
            case 3:
                printf("Using infinite plane wave\n\n");
                break;
            case 4:
                printf("Using pencil beam\n\n");
                break;
            case 5:
                printf("Using top-hat focus, Gaussian far field beam\n\n");
                break;
            case 6:
                printf("Using Gaussian focus, top-hat far field beam\n\n");
                break;
            case 7:
                printf("Using Laguerre-Gaussian LG01 focus, LG01 far field beam\n\n");
                break;
        }

        printf("xFocus = %0.8f [cm]\n",xFocus);
        printf("yFocus = %0.8f [cm]\n",yFocus);
        printf("zFocus = %0.8f [cm]\n",zFocus);

        if (beamtypeFlag!=2) {
            printf("Beam direction defined by unit vector:\n");
            printf("ux0 = %0.8f\n",ux0);
            printf("uy0 = %0.8f\n",uy0);
            printf("uz0 = %0.8f\n",uz0);
        }
        if (beamtypeFlag<=1) {
            printf("Beam waist radius = %0.8f [cm]\n",waist);
            printf("Beam divergence half-angle = %0.8f [rad]\n",divergence);
        }
    }
    
    if (boundaryFlag==0) printf("boundaryFlag = 0, so no boundaries.\n");
    else if (boundaryFlag==1) printf("boundaryFlag = 1, so escape at all boundaries.\n");    
	else if (boundaryFlag==2) printf("boundaryFlag = 2, so escape at surface only.\n");    
	else {
        printf("Error: Invalid boundaryFlag\n");
        return;
    }
	
    printf("\n%d tissues used:\n\n",nT);
    for (i=0; i<nT; i++) {
        printf("muav[%ld] = %0.8f [cm^-1]\n",i+1,muav[i]); // +1 because we convert from C's 0-based indexing to MATLAB's 1-based indexing
        printf("musv[%ld] = %0.8f [cm^-1]\n",i+1,musv[i]); // +1 because we convert from C's 0-based indexing to MATLAB's 1-based indexing
        printf("  gv[%ld] = %0.8f [--]\n\n",i+1,gv[i]); // +1 because we convert from C's 0-based indexing to MATLAB's 1-based indexing
    }
    
    // Show tissue on screen, along central z-axis, by listing tissue type #'s.
    printf("Central axial profile of tissue types:\n");
    for (i=0; i<nz; i++) {
        printf("%d",v[i*ny*nx + ny/2*nx + nx/2]+1); // +1 because we convert from C's 0-based indexing to MATLAB's 1-based indexing
    }
    printf("\n\n");
    
	/**************************
	 * ============================ MAJOR CYCLE ========================
	 **********/
	clock_gettime(CLOCK_MONOTONIC, &simulationTimeStart);
	currentCalendarTime = time(NULL);
	printf("%s", ctime(&currentCalendarTime));	
	
	// INITIALIZATIONS
	if (sourceDist) {
        S_CDF = mxMalloc((nx*ny*nz + 1)*sizeof(double));
        S_CDF[0] = 0;
        for(i=1;i<(nx*ny*nz+1);i++) S_CDF[i] = S_CDF[i-1] + S[i-1];
        sourcesum = S_CDF[nx*ny*nz];
        for(i=1;i<(nx*ny*nz+1);i++) S_CDF[i] /= sourcesum;
    } else if (uz0!=1) {
		vx0 = -uy0/sqrt(uy0*uy0 + ux0*ux0);
		vy0 = ux0/sqrt(uy0*uy0 + ux0*ux0);
		vz0 = 0;
	}

    printf("Simulation time requested = %0.2f min\n\n",simulationTimeRequested);
	
	/**** RUN
	 Launch photons, initializing each one before propagation.
	 *****/
	printf("------------- Begin Monte Carlo -------------\n");
    mexEvalString("drawnow;");
	
    nPhotons = 0;
	
    #ifdef _WIN32
    #pragma omp parallel num_threads(omp_get_num_procs())
    #endif
	{
		double	x = 0, y = 0, z = 0;        /* photon position */
		double	ux, uy, uz;     /* photon trajectory as unit vector composants */
        double  rand;           /* Random number used for lookup into CDF of 3D source */
        long    d;              /* Used in CDF lookup */
		double  ux_temp;        /* temporary value used during SPIN */
		double	s;              /* step sizes. s = -log(RND)/mus [cm] */
		double  stepLeft;       /* dimensionless step size*/
		double	costheta;       /* cos(theta) */
		double  sintheta;       /* sin(theta) */
		double	cospsi;         /* cos(psi) */
		double  sinpsi;         /* sin(psi) */
		double	psi;            /* azimuthal angle */
		double	photonWeight;   /* photon weight */
		double	absorb;         /* weighted deposited in a step due to absorption */
		bool    photonAlive;    /* flag, true or false */
		bool    sameVoxel;      /* Are they in the same voxel? */
		double	mua;            /* absorption coefficient [cm^-1] */
		double	mus;            /* scattering coefficient [cm^-1] */
		double	g;              /* anisotropy [-] */
		double  xTarget, yTarget, zTarget;
		double	wx0, wy0, wz0;  /* normal vector 2 to photon trajectory. Not necessarily normal to v0. */
		double	r, theta, phi;
        long    j;        /* General-purpose thread-specific index variable */
		double	xTentative, yTentative, zTentative; /* temporary variables, used during photon step. */
		long 	ix, iy, iz;     /* Used to track photons */
		bool    photonInsideVolume;  /* flag, true or false */
		dsfmt_t dsfmt;			/* Thread-specific "state" of dSFMT pseudo-random number generator */

        #ifdef _WIN32
		dsfmt_init_gen_rand(&dsfmt,simulationTimeStart.tv_nsec + omp_get_thread_num()); // Seed random number generator
        //printf("Thread #%i seeded with %li\n",omp_get_thread_num()+1,(uint32_t)(simulationTimeStart.tv_nsec + omp_get_thread_num()));
        #else
		dsfmt_init_gen_rand(&dsfmt,simulationTimeStart.tv_nsec); // Seed random number generator
        #endif
		do {
			/**** LAUNCH 
			 Initialize photon position and trajectory.
			 *****/

			photonWeight = 1;
			photonAlive = true;
			sameVoxel = false; // Photon is initialized as if it has just entered the voxel it's created in, to ensure proper initialization of voxel index and ray properties
			
            #ifdef _WIN32
			#pragma omp atomic
            #endif
            nPhotons += 1;				/* increment photon count */
			
			/****************************/
			/* Initial position and trajectory */
            if (sourceDist) {
                // Search the cumulative distribution function via binary tree method to find the voxel to start the photon in
                rand = RandomNum;
                d = nx*ny*nz-1; // Number of elements in the current section to be searched, minus one
                j = d/2; // Index of the middle element of the current section
                while (!(S_CDF[j] < rand && S_CDF[j+1] >= rand)) { // Binary tree search
                    if (S_CDF[j] >= rand) {
                        j = j + (d-2)/4 - d/2;
                        d = d/2 - 1;
                    } else {
                        d = d - d/2 - 1;
                        j = j + d/2 + 1;
                    }
                }
                
                ix = j%nx;
                iy = j/nx%ny;
                iz = j/nx/ny;
                
                x = (ix - nx/2.0 + 1 - RandomNum)*dx;
                y = (iy - ny/2.0 + 1 - RandomNum)*dy;
                z = (iz          + 1 - RandomNum)*dz;
                
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

                    x		= xTarget - zTarget/uz*ux; // the coordinates for the ray starting point is the intersection of the ray with the z = 0 surface
                    y		= yTarget - zTarget/uz*uy;
                    z		= 0;
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

                    x		= xTarget - zTarget/uz*ux; // the coordinates for the ray starting point is the intersection of the ray with the z = 0 surface
                    y		= yTarget - zTarget/uz*uy;
                    z		= 0;
                break;
                case 2: // isotropically emitting point source
                    costheta = 1 - 2*RandomNum;
                    sintheta = sqrt(1 - costheta*costheta);
                    psi = 2*PI*RandomNum;
                    x = xFocus;
                    y = yFocus;
                    z = zFocus;
                    ux = sintheta*cos(psi);
                    uy = sintheta*sin(psi);
                    uz = costheta;
                break;
                case 3: // infinite plane wave
                    if (boundaryFlag==1) {
                        x = nx*dx*(RandomNum-0.5); // Generates a random x coordinate within the cuboid
                        y = ny*dy*(RandomNum-0.5); // Generates a random y coordinate within the cuboid
                    }
                    else {
                        x = infWaveScaleFactor*nx*dx*(RandomNum-0.5); // Generates a random x coordinate within an interval infWaveScaleFactor times the cuboid's size
                        y = infWaveScaleFactor*ny*dy*(RandomNum-0.5); // Generates a random y coordinate within an interval infWaveScaleFactor times the cuboid's size
                    }
                    z = 0;
                    ux = ux0;
                    uy = uy0;
                    uz = uz0;
                break;
                case 4: // pencil beam
                    x	= xFocus - zFocus/uz0*ux0; 
                    y	= yFocus - zFocus/uz0*uy0;
                    z	= 0;
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

                    x		= xTarget - zTarget/uz*ux; // the coordinates for the ray starting point is the intersection of the ray with the z = 0 surface
                    y		= yTarget - zTarget/uz*uy;
                    z		= 0;
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

                    x		= xTarget - zTarget/uz*ux; // the coordinates for the ray starting point is the intersection of the ray with the z = 0 surface
                    y		= yTarget - zTarget/uz*uy;
                    z		= 0;
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

                    x		= xTarget - zTarget/uz*ux; // the coordinates for the ray starting point is the intersection of the ray with the z = 0 surface
                    y		= yTarget - zTarget/uz*uy;
                    z		= 0;
                break;
            }

			/* HOP_DROP_SPIN_CHECK
			 Propagate one photon until it dies as determined by ROULETTE.
			 *******/
			do { // while photonAlive
				/**** HOP
				 Take step to new position
				 *****/
				stepLeft	= -log(RandomNum);				/* dimensionless step */
				
				do {  // while stepLeft>0
					if (!sameVoxel) {
						ix = floor(nx/2.0 + x/dx);
						iy = floor(ny/2.0 + y/dy);
						iz = floor(z/dz);        
						
						photonInsideVolume = true;  // Initialize the photon as inside the volume, then check.
						if (boundaryFlag==0) { // Infinite medium.
									// Check if photon has wandered outside volume.
									// If so, set tissue type to boundary value, but let photon wander.
									// Set photonInsideVolume to false, so DROP does not deposit energy.
							if (iz>=nz) {iz=nz-1; photonInsideVolume = false;}
							if (ix>=nx) {ix=nx-1; photonInsideVolume = false;}
							if (iy>=ny) {iy=ny-1; photonInsideVolume = false;}
							if (iz<0)   {iz=0;    photonInsideVolume = false;}
							if (ix<0)   {ix=0;    photonInsideVolume = false;}
							if (iy<0)   {iy=0;    photonInsideVolume = false;}
						} else if (boundaryFlag==1) { // Escape at boundaries
							if (iz>=nz) {iz=nz-1; photonAlive = false;}
							if (ix>=nx) {ix=nx-1; photonAlive = false;}
							if (iy>=ny) {iy=ny-1; photonAlive = false;}
							if (iz<0)   {iz=0;    photonAlive = false;}
							if (ix<0)   {ix=0;    photonAlive = false;}
							if (iy<0)   {iy=0;    photonAlive = false;}
						} else if (boundaryFlag==2) { // Escape at top surface, no x,y bottom z boundaries
							if (iz>=nz) {iz=nz-1; photonInsideVolume = false;}
							if (ix>=nx) {ix=nx-1; photonInsideVolume = false;}
							if (iy>=ny) {iy=ny-1; photonInsideVolume = false;}
							if (iz<0)   {iz=0;    photonAlive = false;}
							if (ix<0)   {ix=0;    photonInsideVolume = false;}
							if (iy<0)   {iy=0;    photonInsideVolume = false;}
						}
                        if (!photonInsideVolume) {
                            if ((fabs(x) > infWaveScaleFactor*nx*dx/2.0) || (fabs(y) > infWaveScaleFactor*ny*dy/2.0) || (fabs(z-nz*dz/2.0) > infWaveScaleFactor*nz*dz/2.0)) {photonAlive = false;};
                        }
						if (!photonAlive) break;
						
						/* Get tissue voxel properties of current position.
						 * If photon beyond outer edge of defined voxels, 
						 * the tissue equals properties of outermost voxels.
						 */
						j		= iz*ny*nx + iy*nx + ix;
						mua 	= muav[v[j]];
						mus 	= musv[v[j]];
						g 		= gv[v[j]];
					}
					
					s     = stepLeft/mus;				/* Step size [cm].*/
					xTentative = x + s*ux;				/* Update positions. [cm] */
					yTentative = y + s*uy;	
					zTentative = z + s*uz;
					
                    // Now determine if the tentative position is in the same voxel.
                    // The nx/2.0 and ny/2.0 terms are necessary to take into account the offset of voxel surfaces in the case of odd nx or ny.
                    sameVoxel = floor(nx/2.0 + x/dx) == floor(nx/2.0 + xTentative/dx) &&
                                floor(ny/2.0 + y/dy) == floor(ny/2.0 + yTentative/dy) &&
                                floor(         z/dz) == floor(         zTentative/dz);
                    
					if (sameVoxel) { /* photon in same voxel */
						x=xTentative;					/* Update positions. */
						y=yTentative;
						z=zTentative;
						
						/**** DROP
						 Drop photon weight (photonWeight) into local bin.
						 *****/
						absorb = photonWeight*(1 - exp(-mua*s));   /* photon weight absorbed at this step */
						photonWeight -= absorb;					   /* decrement WEIGHT by amount absorbed */
						if (photonInsideVolume) {	// only save data if the photon is inside simulation cuboid
                            #ifdef _WIN32
							#pragma omp atomic
                            #endif
								F[j] += absorb;
						}
						
						stepLeft = 0;		/* dimensionless step remaining */
                    } else { /* photon has crossed voxel boundary */
						/* step to voxel face + "littlest step" so just inside new voxel. */
						s = ls + FindVoxelFace(x,y,z,dx,dy,dz,nx,ny,ux,uy,uz);
						
						/**** DROP
						 Drop photon weight (photonWeight) into local bin.
						 *****/
						absorb = photonWeight*(1 - exp(-mua*s));   /* photon weight absorbed at this step */
						photonWeight -= absorb;                    /* decrement WEIGHT by amount absorbed */
						if (photonInsideVolume) {	// only save data if the photon is inside simulation cuboid
                            #ifdef _WIN32
							#pragma omp atomic
                            #endif
								F[j] += absorb;
						}
						
						stepLeft -= s*mus;  /* dimensionless step remaining */
						if (stepLeft<=ls) stepLeft = 0;
						
						x += s*ux;
						y += s*uy;
						z += s*uz;
					}
				} while(stepLeft>0);
				
				/**** CHECK ROULETTE 
				 If photon weight below THRESHOLD, then terminate photon using Roulette technique.
				 Photon has CHANCE probability of having its weight increased by factor of 1/CHANCE,
				 and 1-CHANCE probability of terminating.
				 *****/
				if (photonWeight < THRESHOLD) {
					if (RandomNum <= CHANCE) photonWeight /= CHANCE;
					else photonAlive = false;
				}
				
				if (photonAlive) {
					/**** SPIN 
					 Scatter photon into new trajectory defined by theta and psi.
					 Theta is specified by cos(theta), which is determined 
					 based on the Henyey-Greenstein scattering function.
					 Convert theta and psi into cosines ux, uy, uz. 
					 *****/
					/* Sample for costheta */
                    costheta = g? (1 + g*g - pow((1 - g*g)/(1 - g + 2*g*RandomNum),2))/(2*g) : 2*RandomNum - 1;
					sintheta = sqrt(1 - costheta*costheta);
					
					/* Sample psi. */
					psi = 2*PI*RandomNum;
					cospsi = cos(psi);
					sinpsi = sin(psi);
                    
					/* New trajectory. */
					if (1 - fabs(uz) <= ONE_MINUS_COSZERO) {      /* close to perpendicular. */
						ux = sintheta*cospsi;
						uy = sintheta*sinpsi;
						uz = costheta*SIGN(uz);
					} else {					                  /* usually use this option */
						ux_temp = sintheta*(ux*uz*cospsi - uy*sinpsi)/sqrt(1 - uz*uz) + ux*costheta;
						uy      = sintheta*(uy*uz*cospsi + ux*sinpsi)/sqrt(1 - uz*uz) + uy*costheta;
						uz      = -sintheta*cospsi*sqrt(1 - uz*uz) + uz*costheta;
                        ux      = ux_temp;
					}
				}
			} while (photonAlive);
			
            #ifdef _WIN32
			#pragma omp master
            #endif
			{
				// Check whether ctrl+c has been pressed
                #ifdef _WIN32
                if(utIsInterruptPending()) {
                    ctrlc_caught = true;
                    printf("Ctrl+C detected, stopping.\n");
                }
                #endif
                clock_gettime(CLOCK_MONOTONIC, &simulationTimeCurrent);
 
				// Print out message about progress.
                newPctProgress = 100*(simulationTimeCurrent.tv_sec - simulationTimeStart.tv_sec + (simulationTimeCurrent.tv_nsec - simulationTimeStart.tv_nsec)/1000000000.0) / (simulationTimeRequested*60);
				if ((newPctProgress != pctProgress) && !ctrlc_caught) {
					pctProgress = newPctProgress;
					if ((pctProgress<10) || (pctProgress>90) || (fmod(pctProgress, 10)==0)) {
						printf("%i%% done\n", pctProgress);
                        mexEvalString("drawnow;");
					}
				}
			}
		} while ((pctProgress < 100) && !ctrlc_caught);  /* end RUN */
	}
    
	printf("------------------------------------------------------\n");
	clock_gettime(CLOCK_MONOTONIC, &simulationTimeCurrent);
	simulationTimeSpent = simulationTimeCurrent.tv_sec - simulationTimeStart.tv_sec + (simulationTimeCurrent.tv_nsec - simulationTimeStart.tv_nsec)/1000000000.0;
	printf("Simulated %0.1e photons at a rate of %0.1e photons per minute\n",(double)nPhotons, nPhotons*60/simulationTimeSpent);
	
    if (sourceDist) mxFree(S_CDF);
    
    /**** SAVE
     Convert data to relative fluence rate [cm^-2] and return.
     *****/
    
    // Normalize deposition to yield fluence rate (F).
    if (sourceDist) {
		for (i=0; i<nx*ny*nz;i++){
			F[i] /= nPhotons*muav[v[i]]/sourcesum; // Absolute fluence rate [W/cm^2]
		}
    } else if (beamtypeFlag == 3 && boundaryFlag != 1) {
		for (i=0; i<nx*ny*nz;i++){
			F[i] /= dx*dy*dz*nPhotons/(infWaveScaleFactor*infWaveScaleFactor)*muav[v[i]]; // Normalized fluence rate [W/cm^2/W.incident]
		}
	} else {
		for (i=0; i<nx*ny*nz;i++){
			F[i] /= dx*dy*dz*nPhotons*muav[v[i]]; // Normalized fluence rate [W/cm^2/W.incident]
		}
	}
	
    printf("------------------------------------------------------\n");
    currentCalendarTime = time(NULL);
    printf("%s\n", ctime(&currentCalendarTime));

    return;
} /* end of main */
