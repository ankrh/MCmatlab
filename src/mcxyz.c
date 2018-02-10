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
 *
 ** COMPILING ON WINDOWS
 * Can be compiled in MATLAB with "mex COPTIMFLAGS='$COPTIMFLAGS -Ofast -fopenmp' LDOPTIMFLAGS='$LDOPTIMFLAGS -Ofast -fopenmp' .\src\mcxyz.c ".\src\libut.lib""
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
 * Compile in MATLAB with "mex COPTIMFLAGS='$COPTIMFLAGS -Ofast' LDOPTIMFLAGS='$LDOPTIMFLAGS -Ofast' ./src/mcxyz.c"
 *
 * To get the MATLAB C compiler to work, try this:
 * 1. Install XCode from the App Store
 * 2. Type "mex -setup" in the MATLAB command window
 ********************************************/

#include "mex.h"
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdbool.h>
#include <time.h>
#ifdef _WIN32 // This is defined on both win32 and win64 systems, we use this preprocessor condition to avoid loading openmp or libut on, e.g., Mac
#include <omp.h>
#endif
#define DSFMT_MEXP 19937 // Mersenne exponent for dSFMT
#include "dSFMT-src-2.2.3/dSFMT.c" //  double precision SIMD oriented Fast Mersenne Twister(dSFMT)

#define Ntiss		19          /* Number of tissue types. */
#define STRLEN 		32          /* String length. */
#define ls          1.0E-7      /* Moving photon a little bit off the voxel face */
#define	PI          3.1415926
#define THRESHOLD   0.01		/* used in roulette */
#define CHANCE      0.1  		/* used in roulette */
#define SQR(x)		(x*x) 
#define SIGN(x)     ((x)>=0 ? 1:-1)
#define RandomNum   dsfmt_genrand_open_close(&dsfmt) /* Calls for a random number in (0,1] */
#define ONE_MINUS_COSZERO 1.0E-12   /* If 1-cos(theta) <= ONE_MINUS_COSZERO, fabs(theta) <= 1e-6 rad. */
/* If 1+cos(theta) <= ONE_MINUS_COSZERO, fabs(PI-theta) <= 1e-6 rad. */

/* DECLARE FUNCTIONS */
bool SameVoxel(double x1,double y1,double z1, double x2, double y2, double z2, double dx,double dy,double dz);
/* Asks,"In the same voxel?" */
double min3(double a, double b, double c);
double FindVoxelFace(double x1,double y1,double z1, double x2, double y2, double z2,double dx,double dy,double dz, double ux, double uy, double uz);
/* How much step size will the photon take to get the first voxel crossing in one single long step? */
void axisrotate(double* x,double* y,double* z, double ux, double uy, double uz,double theta);
/* Rotates a point with coordinates (x,y,z) an angle theta around axis with direction vector (ux,uy,uz) */

void mexFunction( int nlhs, mxArray *plhs[],
        int nrhs, mxArray const *prhs[] ) {
    
#ifdef _WIN32
extern bool utIsInterruptPending(); // Allows catching ctrl+c while executing the mex function
#endif

    // Extracting quantities from MCinput struct
    double	simulationTimeRequested     =      *mxGetPr(mxGetField(prhs[0],0,"simulationTime"));
    int		beamtypeFlag                = (int)*mxGetPr(mxGetField(prhs[0],0,"beamtypeFlag"));
    int		boundaryFlag                = (int)*mxGetPr(mxGetField(prhs[0],0,"boundaryFlag"));
    double	xFocus                      =      *mxGetPr(mxGetField(prhs[0],0,"xFocus"));
    double	yFocus                      =      *mxGetPr(mxGetField(prhs[0],0,"yFocus"));
    double	zFocus                      =      *mxGetPr(mxGetField(prhs[0],0,"zFocus"));
    double	ux0                         =      *mxGetPr(mxGetField(prhs[0],0,"ux0"));
    double	uy0                         =      *mxGetPr(mxGetField(prhs[0],0,"uy0"));
    double	uz0                         =      *mxGetPr(mxGetField(prhs[0],0,"uz0"));
    double	waist                       =      *mxGetPr(mxGetField(prhs[0],0,"waist"));
    double	divergence                  =      *mxGetPr(mxGetField(prhs[0],0,"divergence"));
    double	dx                          =      *mxGetPr(mxGetField(prhs[0],0,"dx"));
    double	dy                          =      *mxGetPr(mxGetField(prhs[0],0,"dy"));
    double	dz                          =      *mxGetPr(mxGetField(prhs[0],0,"dz"));
    
    mxArray *tissueList = mxGetField(prhs[0],0,"tissueList");
    mxArray *T          = mxGetField(prhs[0],0,"T");
    
    mwSize const *dimPtr = mxGetDimensions(T);
    int     nx = dimPtr[0];
    int     ny = dimPtr[1];
    int     nz = dimPtr[2];
    int     nT = mxGetN(tissueList);
    char    *v = (char *)mxGetData(T);
    
    long	i;
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
	unsigned long long	nPhotons;           /* number of photons in simulation */
	int     pctProgress, newPctProgress;    /* Simulation progress in percent */
    bool    ctrlc_caught = false;           // Has a ctrl+c been passed from MATLAB?
	double  extendedBoxScaleFactor;
	extendedBoxScaleFactor = 6.0;
    /* The extendedBoxScaleFactor determines the size of the region that 
     * photons are launched in for the infinite plane wave without 
     * boundaries but also the region that photons are allowed to stay 
     * alive in (if outside, the probability of returning to the region
     * of interest is judged as too low)
     */
    
	/* launch parameters */
	double	vx0, vy0, vz0;                  /* normal vector to photon trajectory */
	
    /* time */
	double	simulationTimeSpent;               // Requested and spent time durations of computation.
	time_t	currentCalendarTime;
	struct  timespec simulationTimeStart, simulationTimeCurrent;
	
    // Display parameters
    printf("nx = %d, dx = %0.8f [cm]\n",nx,dx);
    printf("ny = %d, dy = %0.8f [cm]\n",ny,dy);
    printf("nz = %d, dz = %0.8f [cm]\n",nz,dz);

    printf("beamtypeFlag = %d, ",beamtypeFlag);
    switch(beamtypeFlag) {
		case 0:
			printf("launching top-hat focus, top-hat far field beam\n");
			break;
		case 1:
			printf("launching Gaussian focus, Gaussian far field beam\n");
			break;
		case 2:
			printf("launching isotropic point source\n");
			break;
		case 3:
			printf("launching infinite plane wave\n");
			break;
		case 4:
			printf("launching pencil beam\n");
			break;
		case 5:
			printf("launching top-hat focus, Gaussian far field beam\n");
			break;
		case 6:
			printf("launching Gaussian focus, top-hat far field beam\n");
			break;
	}

    if (boundaryFlag==0)
		printf("boundaryFlag = 0, so no boundaries.\n");
    else if (boundaryFlag==1)
		printf("boundaryFlag = 1, so escape at all boundaries.\n");    
	else if (boundaryFlag==2)
		printf("boundaryFlag = 2, so escape at surface only.\n");    
	else{
        printf("improper boundaryFlag. quit.\n");
        return;
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
		printf("beam waist radius = %0.8f [cm]\n",waist);
		printf("beam divergence half-angle = %0.8f [rad]\n",divergence);
	}
    printf("\n%d tissues available:\n\n",nT);
    for (i=0; i<nT; i++) {
        printf("muav[%ld] = %0.8f [cm^-1]\n",i+1,muav[i]); // +1 because we convert from C's 0-based indexing to MATLAB's 1-based indexing
        printf("musv[%ld] = %0.8f [cm^-1]\n",i+1,musv[i]); // +1 because we convert from C's 0-based indexing to MATLAB's 1-based indexing
        printf("  gv[%ld] = %0.8f [--]\n\n",i+1,gv[i]); // +1 because we convert from C's 0-based indexing to MATLAB's 1-based indexing
    }
    
    // Show tissue on screen, along central z-axis, by listing tissue type #'s.
    printf("Central axial profile of tissue types:\n");
    for (i=0; i<nz; i++) {
        printf("%d",v[(long)(i*ny*nx + ny/2*nx + nx/2)]+1); // +1 because we convert from C's 0-based indexing to MATLAB's 1-based indexing
    }
    printf("\n\n");
    
	/**************************
	 * ============================ MAJOR CYCLE ========================
	 **********/
	clock_gettime(CLOCK_MONOTONIC, &simulationTimeStart);
	currentCalendarTime = time(NULL);
	printf("%s", ctime(&currentCalendarTime));	
	
	/**** INITIALIZATIONS 
	 *****/
	
	if (uz0==1) {
		vx0 = 1;
		vy0 = 0;
		vz0 = 0;
	}
	else {
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
		double	x, y, z;        /* photon position */
		double	ux, uy, uz;     /* photon trajectory as unit vector composants */
		double  uxx, uyy, uzz;	/* temporary values used during SPIN */
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
		long 	j;
		double	xTentative, yTentative, zTentative; /* temporary variables, used during photon step. */
		int 	ix, iy, iz;     /* Added. Used to track photons */
		bool    photonInsideVolume;  /* flag, true or false */
		//int		photonStepsCounter;
		int 	type;
		dsfmt_t dsfmt;			/* Thread-specific "state" of dSFMT pseudo-random number generator */

        #ifdef _WIN32
		dsfmt_init_gen_rand(&dsfmt,(unsigned long)(simulationTimeStart.tv_nsec + omp_get_thread_num())); // Seed random number generator
        //printf("Thread #%i seeded with %li\n",omp_get_thread_num()+1,(unsigned long)(simulationTimeStart.tv_nsec + omp_get_thread_num()));
        #else
		dsfmt_init_gen_rand(&dsfmt,(unsigned long)(simulationTimeStart.tv_nsec)); // Seed random number generator
        #endif
		do {
			/**** LAUNCH 
			 Initialize photon position and trajectory.
			 *****/

			photonWeight = 1.0;                    /* set photon weight to one */
			photonAlive = true;      /* Launch an ALIVE photon */
			sameVoxel = false;					/* Photon is initialized as if it has just entered the voxel it's created in, to ensure proper initialization of voxel index and ray properties */
			//photonStepsCounter = 0;
			
            #ifdef _WIN32
			#pragma omp atomic
            #endif
				nPhotons += 1;				/* increment photon count */
			
			/****************************/
			/* Initial position and trajectory */
			if (beamtypeFlag==0) { // top-hat focus, top-hat far field beam
				r		= waist*sqrt(RandomNum); // for target calculation
				phi		= RandomNum*2.0*PI;
				wx0		= vx0;
				wy0		= vy0;
				wz0		= vz0;
				axisrotate(&wx0,&wy0,&wz0,ux0,uy0,uz0,phi); // w0 unit vector now points in the direction from focus center point to ray target point
				
				xTarget	= xFocus + r*wx0;
				yTarget	= yFocus + r*wy0;
				zTarget = zFocus + r*wz0;
				
				theta	= divergence*sqrt(RandomNum); // for trajectory calculation. The sqrt is valid within paraxial approximation.
				phi		= RandomNum*2.0*PI;
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
			}
			else if (beamtypeFlag==1) { // Gaussian focus, Gaussian far field beam
				r		= waist*sqrt(-0.5*log(RandomNum)); // for target calculation
				phi		= RandomNum*2.0*PI;
				wx0		= vx0;
				wy0		= vy0;
				wz0		= vz0;
				axisrotate(&wx0,&wy0,&wz0,ux0,uy0,uz0,phi); // w0 unit vector now points in the direction from focus center point to ray target point
				
				xTarget	= xFocus + r*wx0;
				yTarget	= yFocus + r*wy0;
				zTarget = zFocus + r*wz0;
				
				theta	= divergence*sqrt(-0.5*log(RandomNum)); // for trajectory calculation. The sqrt is valid within paraxial approximation.
				phi		= RandomNum*2.0*PI;
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
			}
			else if (beamtypeFlag==2) { // isotropically emitting point source
				costheta = 1.0 - 2.0*RandomNum;
				sintheta = sqrt(1.0 - costheta*costheta);
				psi = 2.0*PI*RandomNum;
				cospsi = cos(psi);
				if (psi < PI)
					sinpsi = sqrt(1.0 - cospsi*cospsi); 
				else
					sinpsi = -sqrt(1.0 - cospsi*cospsi);
				x = xFocus;
				y = yFocus;
				z = zFocus;
				ux = sintheta*cospsi;
				uy = sintheta*sinpsi;
				uz = costheta;
			}
			else if (beamtypeFlag==3) { // infinite plane wave
				if (boundaryFlag==1) {
					x = nx*dx*(RandomNum-0.5); // Generates a random x coordinate within the box
					y = ny*dy*(RandomNum-0.5); // Generates a random y coordinate within the box
				}
				else {
					x = extendedBoxScaleFactor*nx*dx*(RandomNum-0.5); // Generates a random x coordinate within an interval extendedBoxScaleFactor times the box' size
					y = extendedBoxScaleFactor*ny*dy*(RandomNum-0.5); // Generates a random y coordinate within an interval extendedBoxScaleFactor times the box' size
				}
				z = 0;
				ux = ux0;
				uy = uy0;
				uz = uz0;
			}
			else if (beamtypeFlag==4) { // pencil beam
				x	= xFocus - zFocus/uz0*ux0; 
				y	= yFocus - zFocus/uz0*uy0;
				z	= 0;
				ux	= ux0;
				uy	= uy0;
				uz	= uz0;
			}
			else if (beamtypeFlag==5) { // top-hat focus, Gaussian far field beam
				r		= waist*sqrt(RandomNum); // for target calculation
				phi		= RandomNum*2.0*PI;
				wx0		= vx0;
				wy0		= vy0;
				wz0		= vz0;
				axisrotate(&wx0,&wy0,&wz0,ux0,uy0,uz0,phi); // w0 unit vector now points in the direction from focus center point to ray target point
				
				xTarget	= xFocus + r*wx0;
				yTarget	= yFocus + r*wy0;
				zTarget = zFocus + r*wz0;
				
				theta	= divergence*sqrt(-0.5*log(RandomNum)); // for trajectory calculation. The sqrt is valid within paraxial approximation.
				phi		= RandomNum*2.0*PI;
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
			}
			else if (beamtypeFlag==6) { // Gaussian focus, top-hat far field beam
				r		= waist*sqrt(-0.5*log(RandomNum)); // for target calculation
				phi		= RandomNum*2.0*PI;
				wx0		= vx0;
				wy0		= vy0;
				wz0		= vz0;
				axisrotate(&wx0,&wy0,&wz0,ux0,uy0,uz0,phi); // w0 unit vector now points in the direction from focus center point to ray target point
				
				xTarget	= xFocus + r*wx0;
				yTarget	= yFocus + r*wy0;
				zTarget = zFocus + r*wz0;
				
				theta	= divergence*sqrt(RandomNum); // for trajectory calculation. The sqrt is valid within paraxial approximation.
				phi		= RandomNum*2.0*PI;
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
			}

			/****************************/
			
			/* HOP_DROP_SPIN_CHECK
			 Propagate one photon until it dies as determined by ROULETTE.
			 *******/
			do {
				
				/**** HOP
				 Take step to new position
				 s = dimensionless stepsize
				 ux, uy, uz are unit vector components of current photon trajectory
				 *****/
				stepLeft	= -log(RandomNum);				/* dimensionless step */
				//photonStepsCounter += 1;
				
				do{  // while stepLeft>0

					if (!sameVoxel) {
						/* Get tissue voxel properties of current position.
						 * If photon beyond outer edge of defined voxels, 
						 * the tissue equals properties of outermost voxels.
						 * Therefore, set outermost voxels to infinite background value.
						 */
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
						}
						else if (boundaryFlag==1) { // Escape at boundaries
							if (iz>=nz) {iz=nz-1; photonAlive = false;}
							if (ix>=nx) {ix=nx-1; photonAlive = false;}
							if (iy>=ny) {iy=ny-1; photonAlive = false;}
							if (iz<0)   {iz=0;    photonAlive = false;}
							if (ix<0)   {ix=0;    photonAlive = false;}
							if (iy<0)   {iy=0;    photonAlive = false;}
						}
						else if (boundaryFlag==2) { // Escape at top surface, no x,y bottom z boundaries
							if (iz>=nz) {iz=nz-1; photonInsideVolume = false;}
							if (ix>=nx) {ix=nx-1; photonInsideVolume = false;}
							if (iy>=ny) {iy=ny-1; photonInsideVolume = false;}
							if (iz<0)   {iz=0;    photonAlive = false;}
							if (ix<0)   {ix=0;    photonInsideVolume = false;}
							if (iy<0)   {iy=0;    photonInsideVolume = false;}
						}
                        if (!photonInsideVolume) {
                            if (fabs(x) > extendedBoxScaleFactor*nx*dx/2.0 | fabs(y) > extendedBoxScaleFactor*ny*dy/2.0 | fabs(z-nz*dz/2.0) > extendedBoxScaleFactor*nz*dz/2.0) {photonAlive = false;};
                        }
						if (!photonAlive) break;
						
						/* Get the tissue type of located voxel */
						j		= (long)(iz*ny*nx + iy*nx + ix);
						type	= v[j];
						mua 	= muav[type];
						mus 	= musv[type];
						g 		= gv[type];
					}
					
					s     = stepLeft/mus;				/* Step size [cm].*/
					xTentative = x + s*ux;				/* Update positions. [cm] */
					yTentative = y + s*uy;	
					zTentative = z + s*uz;
					
					sameVoxel = SameVoxel(x,y,z, xTentative, yTentative, zTentative, dx,dy,dz);
					if (sameVoxel) /* photon in same voxel */
					{  
						x=xTentative;					/* Update positions. */
						y=yTentative;
						z=zTentative;
						
						/**** DROP
						 Drop photon weight (photonWeight) into local bin.
						 *****/
						absorb = photonWeight*(1 - exp(-mua*s));	/* photon weight absorbed at this step */
						photonWeight -= absorb;					/* decrement WEIGHT by amount absorbed */
						// If photon within volume of heterogeneity, deposit energy in F[]. 
						// Normalize F[] later, when save output. 
						if (photonInsideVolume) {
                            #ifdef _WIN32
							#pragma omp atomic
                            #endif
								F[j] += absorb;	// only save data if the photon is inside simulation cube
						}
						
						/* Update stepLeft */
						stepLeft = 0;		/* dimensionless step remaining */
					}
					else /* photon has crossed voxel boundary */
					{
						/* step to voxel face + "littlest step" so just inside new voxel. */
						s = ls + FindVoxelFace(x,y,z, xTentative,yTentative,zTentative, dx,dy,dz, ux,uy,uz);
						
						/**** DROP
						 Drop photon weight (photonWeight) into local bin.
						 *****/
						absorb = photonWeight*(1-exp(-mua*s));   /* photon weight absorbed at this step */
						photonWeight -= absorb;                  /* decrement WEIGHT by amount absorbed */
						// If photon within volume of heterogeneity, deposit energy in F[]. 
						// Normalize F[] later, when save output. 
						if (photonInsideVolume) {
                            #ifdef _WIN32
							#pragma omp atomic
                            #endif
								F[j] += absorb;	// only save data if the photon is inside simulation cube
						}
						
						/* Update stepLeft */
						stepLeft -= s*mus;  /* dimensionless step remaining */
						if (stepLeft<=ls) stepLeft = 0;
						
						/* Update positions. */
						x += s*ux;
						y += s*uy;
						z += s*uz;
						
					} //(sameVoxel) /* same voxel */
					
				} while(stepLeft>0); //do...while
				
				/**** CHECK ROULETTE 
				 If photon weight below THRESHOLD, then terminate photon using Roulette technique.
				 Photon has CHANCE probability of having its weight increased by factor of 1/CHANCE,
				 and 1-CHANCE probability of terminating.
				 *****/
				if (photonWeight < THRESHOLD) {
					if (RandomNum <= CHANCE)
						photonWeight /= CHANCE;
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
					if (g == 0.0)
						costheta = 2.0*RandomNum - 1.0;
					else {
						costheta = (1.0 + g*g - SQR((1.0 - g*g)/(1.0 - g + 2*g*RandomNum)))/(2.0*g);
					}
					sintheta = sqrt(1.0 - costheta*costheta); /* sqrt() is faster than sin(). */
					
					/* Sample psi. */
					psi = 2.0*PI*RandomNum;
					cospsi = cos(psi);
					if (psi < PI)
						sinpsi = sqrt(1.0 - cospsi*cospsi);     /* sqrt() is faster than sin(). */
					else
						sinpsi = -sqrt(1.0 - cospsi*cospsi);
					
					/* New trajectory. */
					if (1 - fabs(uz) <= ONE_MINUS_COSZERO) {      /* close to perpendicular. */
						uxx = sintheta * cospsi;
						uyy = sintheta * sinpsi;
						uzz = costheta * SIGN(uz);   /* SIGN() is faster than division. */
					} 
					else {					/* usually use this option */
						uxx = sintheta * (ux * uz * cospsi - uy * sinpsi) / sqrt(1.0 - uz * uz) + ux * costheta;
						uyy = sintheta * (uy * uz * cospsi + ux * sinpsi) / sqrt(1.0 - uz * uz) + uy * costheta;
						uzz = -sintheta * cospsi * sqrt(1.0 - uz * uz) + uz * costheta;
					}
					
					/* Update trajectory */
					ux = uxx;
					uy = uyy;
					uz = uzz;
				}
				
			} while (photonAlive);
			/* end STEP_CHECK_HOP_SPIN */
			/* if ALIVE, continue propagating */
			/* If photon DEAD, then launch new photon. */	
			
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
                newPctProgress = 100*(simulationTimeCurrent.tv_sec - simulationTimeStart.tv_sec + (simulationTimeCurrent.tv_nsec - simulationTimeStart.tv_nsec)/1000000000.0) / (simulationTimeRequested*60.0);
				if (newPctProgress != pctProgress & !ctrlc_caught) {
					pctProgress = newPctProgress;
					if ((pctProgress<10) | (pctProgress>90) | (fmod(pctProgress, 10)==0)) {
						printf("%i%% done\n", pctProgress);
                        mexEvalString("drawnow;");
					}
				}
			}
		} while (pctProgress < 100 & !ctrlc_caught);  /* end RUN */
	}
    
	printf("------------------------------------------------------\n");
	clock_gettime(CLOCK_MONOTONIC, &simulationTimeCurrent);
	simulationTimeSpent = (double)(simulationTimeCurrent.tv_sec  - simulationTimeStart.tv_sec ) + (simulationTimeCurrent.tv_nsec - simulationTimeStart.tv_nsec)/1000000000.0;
	printf("Simulated %0.1e photons at a rate of %0.1e photons per minute\n",(double)nPhotons, nPhotons*60/simulationTimeSpent);
	
    /**** SAVE
     Convert data to relative fluence rate [cm^-2] and return.
     *****/
    
    // Normalize deposition to yield fluence rate (F).
	if (beamtypeFlag == 3 && boundaryFlag != 1) {
		for (i=0; i<nx*ny*nz;i++){
			F[i] /= (dx*dy*dz*nPhotons/(extendedBoxScaleFactor*extendedBoxScaleFactor)*muav[v[i]]);
		}
	}
	else {
		for (i=0; i<nx*ny*nz;i++){
			F[i] /= (dx*dy*dz*nPhotons*muav[v[i]]);
		}
	}
	
    printf("------------------------------------------------------\n");
    currentCalendarTime = time(NULL);
    printf("%s\n", ctime(&currentCalendarTime));

    return;
} /* end of main */



/* SUBROUTINES */

/***********************************************************
 *  Determine if the two position are located in the same voxel
 *	Returns 1 if same voxel, 0 if not same voxel.
 ****/				
bool SameVoxel(double x1,double y1,double z1, double x2, double y2, double z2, double dx,double dy,double dz)
{
    double xmin=fmin((floor)(x1/dx),(floor)(x2/dx))*dx;
    double ymin=fmin((floor)(y1/dy),(floor)(y2/dy))*dy;
    double zmin=fmin((floor)(z1/dz),(floor)(z2/dz))*dz;
    double xmax = xmin+dx;
    double ymax = ymin+dy;
    double zmax = zmin+dz;
    bool sameVoxel = false;
    
    sameVoxel=(x1<=xmax && x2<=xmax && y1<=ymax && y2<=ymax && z1<zmax && z2<=zmax);
    return (sameVoxel);
}

/***********************************************************
 * min3
 ****/
double min3(double a, double b, double c) {
    double m;
    if (a <=  fmin(b, c))
        m = a;
    else if (b <= fmin(a, c))
        m = b;
    else
        m = c;
    return m;
}

/********************
 * Steve's version of FindVoxelFace for no scattering.
 * s = ls + FindVoxelFace(x,y,z, xTentative, yTentative, zTentative, dx, dy, dz, ux, uy, uz);
 ****/
double FindVoxelFace(double x1,double y1,double z1, double x2, double y2, double z2,double dx,double dy,double dz, double ux, double uy, double uz)
{	
    int ix1 = floor(x1/dx);
    int iy1 = floor(y1/dy);
    int iz1 = floor(z1/dz);
    
    int ix2,iy2,iz2;
    if (ux>=0)
        ix2=ix1+1;
    else
        ix2 = ix1;
    
    if (uy>=0)
        iy2=iy1+1;
    else
        iy2 = iy1;
    
    if (uz>=0)
        iz2=iz1+1;
    else
        iz2 = iz1;
    
    double xs = fabs( (ix2*dx - x1)/ux);
    double ys = fabs( (iy2*dy - y1)/uy);
    double zs = fabs( (iz2*dz - z1)/uz);
    
    double s = min3(xs,ys,zs);
    
    return (s);
}

/********************
 * Function for rotating a point at coordinates (x,y,z) an angle theta around an axis with direction unit vector (ux,uy,uz).
 ****/
void axisrotate(double* x,double* y,double* z, double ux, double uy, double uz,double theta)
{	
    //printf("Before rotation (x,y,z) = (%0.8f,%0.8f,%0.8f)\n",*x,*y,*z);
    double xin = *x;
	double yin = *y;
	double zin = *z;
	double st = sin(theta);
	double ct = cos(theta);
	
	*x = (ct + ux*ux*(1-ct))*xin		+		(ux*uy*(1-ct) - uz*st)*yin		+		(ux*uz*(1-ct) + uy*st)*zin;
	*y = (uy*ux*(1-ct) + uz*st)*xin		+		(ct + uy*uy*(1-ct))*yin			+		(uy*uz*(1-ct) - ux*st)*zin;
	*z = (uz*ux*(1-ct) - uy*st)*xin		+		(uz*uy*(1-ct) + ux*st)*yin		+		(ct + uz*uz*(1-ct))*zin;
    //printf("After rotation (x,y,z) = (%0.8f,%0.8f,%0.8f)\n",*x,*y,*z);
}
