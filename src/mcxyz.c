/********************************************
 *  mcxyz.c,	in ANSI Standard C programing language
 *      Usage:  mcxyz myname and myname_T.bin
 *      which loads myname_H.mci, and saves myname_F.bin. 
 *
 *	created 2010, 2012 by
 *	Steven L. JACQUES
 *  Ting LI
 *	Oregon Health & Science University
 *
 *  USAGE   mcxyz myname
 *              where myname is the user's choice. 
 *          The program reads two files prepared by user:
 *                  myname_H.mci    = header input file for mcxyz
 *                  myname_T.bin    = tissue structure file
 *          The output will be written to 3 files:
 *                  myname_OP.m     = optical properties  (mua, mus, g for each tissue type)
 *                  myname_F.bin    = fluence rate output F[i] [W/cm^2 per W delivered]
 *
 *  The MATLAB program maketissue.m can create the two input files (myname_H.mci, myname_T.bin).
 *
 *  The MATLAB program lookmcxyz.m can read the output files and display
 *          1. Fluence rate F [W/cm^2 per W delivered]
 *          2. Deposition rate A [W/cm^3 per W delivered].
 *
 *  Log:
 *  Written by Ting based on Steve's mcsub.c., 2010.
 *      Use Ting's FindVoxelFace().
 *	Use Steve's FindVoxelFace(), Dec. 30, 2010.
 *  Reorganized by Steve. May 8, 2012:
 *      Reads input files, outputs binary files.
 *	Edited to included a uniformly distributed light source over the entire surface by Mathias Christensen 09/01/2014
 *  Overhauled by Anders Hansen and Dominik Marti in April 2017. Fundamental method remained unchanged.
 **********/

#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdbool.h>
#include <time.h>
#include <omp.h>
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
double FindVoxelFace2(double x1,double y1,double z1, double x2, double y2, double z2,double dx,double dy,double dz, double ux, double uy, double uz);
/* How much step size will the photon take to get the first voxel crossing in one single long step? */
void axisrotate(double* x,double* y,double* z, double ux, double uy, double uz,double theta);
/* Rotates a point with coordinates (x,y,z) an angle theta around axis with direction vector (ux,uy,uz) */

int main(int argc, const char * argv[]) {
    
    if (argc==0) {
        printf("assuming you've compiled mcxyz.c as mcxyz ...\n");
		printf("USAGE: mcxyz name\n");
		printf("which will load the files name_H.mci and name_T.bin\n");
		printf("and run the Monte Carlo program.\n");
		printf("Yields name_F.bin, which holds the fluence rate distribution.\n");
        return 0;
    }
	
	

	/* other variables */
	unsigned long long	nPhotons;       /* number of photons in simulation */
	int     pctProgress, newPctProgress;    /* Simulation progress in percent */
	long	i,totalVoxelCount;

	/* launch parameters */
	int		beamtypeFlag, boundaryFlag;
	double	xFocus, yFocus, zFocus;
	double	ux0, uy0, uz0;	/* beam center axis direction unit vector composants */
	double	vx0, vy0, vz0;  /* normal vector 1 to photon trajectory */
	double	waist;
	double	divergence;
	double  infPlaneWaveScaleFactor;
	infPlaneWaveScaleFactor = 6.0;
	
	/* mcxyz bin variables */
	double	dx, dy, dz;     /* bin size [cm] */
	int		nx, ny, nz, numberOfTissuetypes; /* # of bins */
    
    /* time */
	double	simulationTimeRequested, simulationTimeSpent;               // Requested and spent time durations of computation.
	time_t	currentCalendarTime;
	struct  timespec simulationTimeStart, simulationTimeCurrent;
	
	/* tissue parameters */
	char	tissuename[50][32];
	double 	muav[Ntiss];            // muav[0:Ntiss-1], absorption coefficient of ith tissue type
	double 	musv[Ntiss];            // scattering coeff. 
	double 	gv[Ntiss];              // anisotropy of scattering
    
	/* Input/Output */
	char   	myname[STRLEN];		// Holds the user's choice of myname, used in input and output files. 
	char	filename[STRLEN];     // temporary filename for writing output.
    FILE*	fid=NULL;               // file ID pointer 
    char    buf[32];                // buffer for reading header.dat
    
    strcpy(myname, argv[1]);    // acquire name from argument of function call by user.
    printf("name = %s\n",myname);
    
	/**** INPUT FILES *****/
    /* IMPORT myname_H.mci */
    strcpy(filename,myname);
    strcat(filename, "_H.mci");
	fid = fopen(filename,"r");
	fgets(buf, 32, fid);
		// run parameters
		sscanf(buf, "%lf", &simulationTimeRequested); // desired time duration of run [min]
		fgets(buf, 32, fid);
		sscanf(buf, "%d", &nx);  // # of bins  
		fgets(buf, 32,fid);
		sscanf(buf, "%d", &ny);  // # of bins
		fgets(buf, 32,fid);
		sscanf(buf, "%d", &nz);  // # of bins   
	
		fgets(buf, 32,fid);
		sscanf(buf, "%lf", &dx);	 // size of bins [cm]
		fgets(buf, 32,fid);
		sscanf(buf, "%lf", &dy);	 // size of bins [cm] 
		fgets(buf, 32,fid);
		sscanf(buf, "%lf", &dz);	 // size of bins [cm] 
	
		// launch parameters
		/*
		beam type: 0 = top-hat focus, top-hat far field beam,
		1 = Gaussian focus, Gaussian far field beam,
		2 = isotropically emitting point, 3 = infinite
		plane wave, 4 = pencil beam, 5 = top-hat focus,
		Gaussian far field beam, 6 = Gaussian focus,
		top-hat far field beam
		*/
		fgets(buf, 32,fid);
		sscanf(buf, "%d", &beamtypeFlag);
		/*
		0 = no boundaries, 1 = escape at boundaries
        2 = escape at surface only. No x, y, bottom z
        boundaries
		*/
        fgets(buf, 32,fid);
        sscanf(buf, "%d", &boundaryFlag);
	
		fgets(buf, 32,fid);
		sscanf(buf, "%lf", &xFocus);  // xFocus
		fgets(buf, 32,fid);
		sscanf(buf, "%lf", &yFocus);  // yFocus
		fgets(buf, 32,fid);
		sscanf(buf, "%lf", &zFocus);  // zFocus

		fgets(buf, 32,fid);
		sscanf(buf, "%lf", &ux0);  // ux trajectory
		fgets(buf, 32,fid);
		sscanf(buf, "%lf", &uy0);  // uy trajectory
		fgets(buf, 32,fid);
		sscanf(buf, "%lf", &uz0);  // uz trajectory

		fgets(buf, 32,fid);
		sscanf(buf, "%lf", &waist);  // waist 1/e^2 radius
		fgets(buf, 32,fid);
		sscanf(buf, "%lf", &divergence);  // divergence 1/e^2 half-angle of incoming beam in rad
    
	
		// tissue optical properties
		fgets(buf, 32,fid);
		sscanf(buf, "%d", &numberOfTissuetypes);				// # of tissue types in tissue list
		for (i=1; i<=numberOfTissuetypes; i++) {
			fgets(buf, 32, fid);
			sscanf(buf, "%lf", &muav[i]);	// absorption coeff [cm^-1]
			fgets(buf, 32, fid);
			sscanf(buf, "%lf", &musv[i]);	// scattering coeff [cm^-1]
			fgets(buf, 32, fid);
			sscanf(buf, "%lf", &gv[i]);		// anisotropy of scatter [dimensionless]
		}    
    fclose(fid);
    
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
        return 0;
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
    printf("\n%d tissues available:\n\n",numberOfTissuetypes);
    for (i=1; i<=numberOfTissuetypes; i++) {
        printf("muav[%ld] = %0.8f [cm^-1]\n",i,muav[i]);
        printf("musv[%ld] = %0.8f [cm^-1]\n",i,musv[i]);
        printf("  gv[%ld] = %0.8f [--]\n\n",i,gv[i]);
    }

    /* IMPORT BINARY TISSUE FILE */
	char 	*v=NULL;
	float 	*F=NULL;
	//float	*R=NULL;
	totalVoxelCount = nx*ny*nz;
	v  = ( char *)malloc(totalVoxelCount*sizeof(char));  /* tissue structure */
	F  = (float *)malloc(totalVoxelCount*sizeof(float));	/* relative fluence rate [W/cm^2/W.delivered] */
	// NOT READY: R  = (float *)malloc(totalVoxelCount*sizeof(float));	/* escaping flux [W/cm^2/W.delivered] */
    
	// read binary file
    strcpy(filename,myname);
    strcat(filename, "_T.bin");
    fid = fopen(filename, "rb");
    fread(v, sizeof(char), totalVoxelCount, fid);
    fclose(fid);
    
    // Show tissue on screen, along central z-axis, by listing tissue type #'s.
    printf("Central axial profile of tissue types:\n");
    for (i=0; i<nz; i++) {
        printf("%d",v[(long)(i*ny*nx + ny/2*nx + nx/2)]);
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
	for(i=0; i<totalVoxelCount;i++) 	F[i] = 0; // ensure F[] starts empty.	
	
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

	nPhotons = 0;
	
	#pragma omp parallel num_threads(omp_get_num_procs())
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

		dsfmt_init_gen_rand(&dsfmt,(unsigned long)(simulationTimeStart.tv_nsec + omp_get_thread_num())); // Seed random number generator
		//printf("Thread #%i seeded with %li\n",omp_get_thread_num()+1,(unsigned long)(simulationTimeStart.tv_nsec + omp_get_thread_num()));
		do {
			/**** LAUNCH 
			 Initialize photon position and trajectory.
			 *****/

			photonWeight = 1.0;                    /* set photon weight to one */
			photonAlive = true;      /* Launch an ALIVE photon */
			sameVoxel = false;					/* Photon is initialized as if it has just entered the voxel it's created in, to ensure proper initialization of voxel index and ray properties */
			//photonStepsCounter = 0;
			
			#pragma omp atomic
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
					x = infPlaneWaveScaleFactor*nx*dx*(RandomNum-0.5); // Generates a random x coordinate within an interval infPlaneWaveScaleFactor times the box' size
					y = infPlaneWaveScaleFactor*ny*dy*(RandomNum-0.5); // Generates a random y coordinate within an interval infPlaneWaveScaleFactor times the box' size
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
							if (iz>=nz) {iz=nz-1; photonAlive = false; stepLeft = 0;}
							if (ix>=nx) {ix=nx-1; photonAlive = false; stepLeft = 0;}
							if (iy>=ny) {iy=ny-1; photonAlive = false; stepLeft = 0;}
							if (iz<0)   {iz=0;    photonAlive = false; stepLeft = 0;}
							if (ix<0)   {ix=0;    photonAlive = false; stepLeft = 0;}
							if (iy<0)   {iy=0;    photonAlive = false; stepLeft = 0;}
						}
						else if (boundaryFlag==2) { // Escape at top surface, no x,y bottom z boundaries
							if (iz>=nz) {iz=nz-1; photonInsideVolume = false;}
							if (ix>=nx) {ix=nx-1; photonInsideVolume = false;}
							if (iy>=ny) {iy=ny-1; photonInsideVolume = false;}
							if (iz<0)   {iz=0;    photonAlive = false; stepLeft = 0;}
							if (ix<0)   {ix=0;    photonInsideVolume = false;}
							if (iy<0)   {iy=0;    photonInsideVolume = false;}
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
							#pragma omp atomic
								F[j] += absorb;	// only save data if the photon is inside simulation cube
						}
						
						/* Update stepLeft */
						stepLeft = 0;		/* dimensionless step remaining */
					}
					else /* photon has crossed voxel boundary */
					{
						/* step to voxel face + "littlest step" so just inside new voxel. */
						s = ls + FindVoxelFace2(x,y,z, xTentative,yTentative,zTentative, dx,dy,dz, ux,uy,uz);
						
						/**** DROP
						 Drop photon weight (photonWeight) into local bin.
						 *****/
						absorb = photonWeight*(1-exp(-mua*s));   /* photon weight absorbed at this step */
						photonWeight -= absorb;                  /* decrement WEIGHT by amount absorbed */
						// If photon within volume of heterogeneity, deposit energy in F[]. 
						// Normalize F[] later, when save output. 
						if (photonInsideVolume) {
							#pragma omp atomic
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
			
			#pragma omp master
			{
				clock_gettime(CLOCK_MONOTONIC, &simulationTimeCurrent);

				// Print out message about progress.
				newPctProgress = 100*(simulationTimeCurrent.tv_sec - simulationTimeStart.tv_sec + (simulationTimeCurrent.tv_nsec - simulationTimeStart.tv_nsec)/1000000000.0) / (simulationTimeRequested*60.0);
				if (newPctProgress != pctProgress) {
					pctProgress = newPctProgress;
					if ((pctProgress<10) | (pctProgress>90) | (fmod(pctProgress, 10)==0)) {
						printf("%i%% done\n", pctProgress);
					}
				}
			}
		} while (pctProgress < 100);  /* end RUN */
	}
    
	printf("------------------------------------------------------\n");
	clock_gettime(CLOCK_MONOTONIC, &simulationTimeCurrent);
	simulationTimeSpent = (double)(simulationTimeCurrent.tv_sec  - simulationTimeStart.tv_sec ) + (simulationTimeCurrent.tv_nsec - simulationTimeStart.tv_nsec)/1000000000.0;
	printf("Simulated %0.1e photons at a rate of %0.1e photons per minute\n",(double)nPhotons, nPhotons*60/simulationTimeSpent);
	
    /**** SAVE
     Convert data to relative fluence rate [cm^-2] and save.
     *****/
    
    // Normalize deposition to yield fluence rate (F).
	if (beamtypeFlag == 3 && boundaryFlag != 1) {
		for (i=0; i<totalVoxelCount;i++){
			F[i] /= (dx*dy*dz*nPhotons/(infPlaneWaveScaleFactor*infPlaneWaveScaleFactor)*muav[v[i]]);
		}
	}
	else {
		for (i=0; i<totalVoxelCount;i++){
			F[i] /= (dx*dy*dz*nPhotons*muav[v[i]]);
		}
	}
    // Save the binary file
    strcpy(filename,myname);
    strcat(filename,"_F.bin");
    printf("Saving %s\n",filename);
    fid = fopen(filename, "wb");   /* 3D voxel output */
    fwrite(F, sizeof(float), totalVoxelCount, fid);
    fclose(fid);
    
    /* save reflectance */
// NOT READY: 
//	strcpy(filename,myname);
//    strcat(filename,"_Ryx.bin");
//    printf("saving %s\n",filename);
//    fid = fopen(filename, "wb");   /* 2D voxel output */
//	int Nyx = ny*nx;
//    fwrite(R, sizeof(float), Nyx, fid);
//    fclose(fid);
//    printf("%s is done.\n",myname);
	
    printf("------------------------------------------------------\n");
    currentCalendarTime = time(NULL);
    printf("%s\n", ctime(&currentCalendarTime));
    
    
    free(v);
    free(F);
    //free(R);
    return 0;
} /* end of main */



/* SUBROUTINES */

/**************************************************************************
 *	RandomGen
 *      A random number generator that generates uniformly
 *      distributed random numbers between 0 and 1 inclusive.
 *      The algorithm is based on:
 *      W.H. Press, S.A. Teukolsky, W.T. Vetterling, and B.P.
 *      Flannery, "Numerical Recipes in C," Cambridge University
 *      Press, 2nd edition, (1992).
 *      and
 *      D.E. Knuth, "Seminumerical Algorithms," 2nd edition, vol. 2
 *      of "The Art of Computer Programming", Addison-Wesley, (1981).
 *
 *      When Type is 0, sets Seed as the seed. Make sure 0<Seed<32000.
 *      When Type is 1, returns a random number.
 *      When Type is 2, gets the status of the generator.
 *      When Type is 3, restores the status of the generator.
 *
 *      The status of the generator is represented by Status[0..56].
 *
 *      Make sure you initialize the seed before you get random
 *      numbers.
 ****/
 
 
 // This is Knuth's subtractive pseudorandom number generator from 1981
#define MBIG 1000000000
#define MSEED 161803398
#define MZ 0
#define FAC 1.0E-9

double RandomGen(char Type, long Seed, long *Status){
	static __thread long i1, i2, ma[56];   /* ma[0] is not used. */
	long        mj, mk;
	short       i, ii;
	
	if (Type == 0) {              /* set seed. */
		//printf("Thread #%i random number generator seed is %ld\n",omp_get_thread_num()+1,Seed);
		mj = MSEED - (Seed < 0 ? -Seed : Seed);
		mj %= MBIG;
		ma[55] = mj;
		mk = 1;
		for (i = 1; i <= 54; i++) {
			ii = (21 * i) % 55;
			ma[ii] = mk;
			mk = mj - mk;
			if (mk < MZ)
				mk += MBIG;
			mj = ma[ii];
		}
		for (ii = 1; ii <= 4; ii++)
			for (i = 1; i <= 55; i++) {
				ma[i] -= ma[1 + (i + 30) % 55];
				if (ma[i] < MZ)
					ma[i] += MBIG;
			}
		i1 = 0;
		i2 = 31;
	} else if (Type == 1) {       /* get a number. */
		if (++i1 == 56)
			i1 = 1;
		if (++i2 == 56)
			i2 = 1;
		mj = ma[i1] - ma[i2];
		if (mj < MZ)
			mj += MBIG;
		ma[i1] = mj;
		return (mj * FAC);
	} else if (Type == 2) {       /* get status. */
		for (i = 0; i < 55; i++)
			Status[i] = ma[i + 1];
		Status[55] = i1;
		Status[56] = i2;
	} else if (Type == 3) {       /* restore status. */
		for (i = 0; i < 55; i++)
			ma[i + 1] = Status[i];
		i1 = Status[55];
		i2 = Status[56];
	} else
		puts("Wrong parameter to RandomGen().");
	return (0);
}
#undef MBIG
#undef MSEED
#undef MZ
#undef FAC


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
 * my version of FindVoxelFace for no scattering.
 * s = ls + FindVoxelFace2(x,y,z, xTentative, yTentative, zTentative, dx, dy, dz, ux, uy, uz);
 ****/
double FindVoxelFace2(double x1,double y1,double z1, double x2, double y2, double z2,double dx,double dy,double dz, double ux, double uy, double uz)
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


/***********************************************************
 *	FRESNEL REFLECTANCE
 * Computes reflectance as photon passes from medium 1 to 
 * medium 2 with refractive indices n1,n2. Incident
 * angle a1 is specified by cosine value ca1 = cos(a1).
 * Program returns value of transmitted angle a1 as
 * value in *ca2_Ptr = cos(a2).
 ****/
double RFresnel(double n1,		/* incident refractive index.*/
                double n2,		/* transmit refractive index.*/
                double ca1,		/* cosine of the incident */
                /* angle a1, 0<a1<90 degrees. */
                double *ca2_Ptr) 	/* pointer to the cosine */
/* of the transmission */
/* angle a2, a2>0. */
{
    double r;
    
    if(n1==n2) { /** matched boundary. **/
        *ca2_Ptr = ca1;
        r = 0.0;
	}
    else if(ca1>(1.0 - 1.0e-12)) { /** normal incidence. **/
        *ca2_Ptr = ca1;
        r = (n2-n1)/(n2+n1);
        r *= r;
	}
    else if(ca1< 1.0e-6)  {	/** very slanted. **/
        *ca2_Ptr = 0.0;
        r = 1.0;
	}
    else  {			  		/** general. **/
        double sa1, sa2; /* sine of incident and transmission angles. */
        double ca2;      /* cosine of transmission angle. */
        sa1 = sqrt(1-ca1*ca1);
        sa2 = n1*sa1/n2;
        if(sa2>=1.0) {	
            /* double check for total internal reflection. */
            *ca2_Ptr = 0.0;
            r = 1.0;
		}
        else {
            double cap, cam;	/* cosines of sum ap or diff am of the two */
            /* angles: ap = a1 + a2, am = a1 - a2. */
            double sap, sam;	/* sines. */
            *ca2_Ptr = ca2 = sqrt(1-sa2*sa2);
            cap = ca1*ca2 - sa1*sa2; /* c+ = cc - ss. */
            cam = ca1*ca2 + sa1*sa2; /* c- = cc + ss. */
            sap = sa1*ca2 + ca1*sa2; /* s+ = sc + cs. */
            sam = sa1*ca2 - ca1*sa2; /* s- = sc - cs. */
            r = 0.5*sam*sam*(cam*cam+cap*cap)/(sap*sap*cam*cam); 
            /* rearranged for speed. */
		}
	}
    return(r);
} /******** END SUBROUTINE **********/



/***********************************************************
 * the boundary is the face of some voxel
 * find the the photon's hitting position on the nearest face of the voxel and update the step size.
 ****/
double FindVoxelFace(double x1,double y1,double z1, double x2, double y2, double z2,double dx,double dy,double dz, double ux, double uy, double uz)
{
    double x_1 = x1/dx;
    double y_1 = y1/dy;
    double z_1 = z1/dz;
    double x_2 = x2/dx;
    double y_2 = y2/dy;
    double z_2 = z2/dz;
    double fx_1 = floor(x_1) ;
    double fy_1 = floor(y_1) ;
    double fz_1 = floor(z_1) ;
    double fx_2 = floor(x_2) ;
    double fy_2 = floor(y_2) ;
    double fz_2 = floor(z_2) ;
    double x=0, y=0, z=0, x0=0, y0=0, z0=0, s=0;
    
    if ((fx_1 != fx_2) && (fy_1 != fy_2) && (fz_1 != fz_2) ) { //#10
        fx_2=fx_1+SIGN(fx_2-fx_1);//added
        fy_2=fy_1+SIGN(fy_2-fy_1);//added
        fz_2=fz_1+SIGN(fz_2-fz_1);//added
        
        x = (fmax(fx_1,fx_2)-x_1)/ux;
        y = (fmax(fy_1,fy_2)-y_1)/uy;
        z = (fmax(fz_1,fz_2)-z_1)/uz;
        if (x == min3(x,y,z)) {
            x0 = fmax(fx_1,fx_2);
            y0 = (x0-x_1)/ux*uy+y_1;
            z0 = (x0-x_1)/ux*uz+z_1;
        }
        else if (y == min3(x,y,z)) {
            y0 = fmax(fy_1,fy_2);
            x0 = (y0-y_1)/uy*ux+x_1;
            z0 = (y0-y_1)/uy*uz+z_1;
        }
        else {
            z0 = fmax(fz_1,fz_2);
            y0 = (z0-z_1)/uz*uy+y_1;
            x0 = (z0-z_1)/uz*ux+x_1;
        }
    }
    else if ( (fx_1 != fx_2) && (fy_1 != fy_2) ) { //#2
        fx_2=fx_1+SIGN(fx_2-fx_1);//added
        fy_2=fy_1+SIGN(fy_2-fy_1);//added
        x = (fmax(fx_1,fx_2)-x_1)/ux;
        y = (fmax(fy_1,fy_2)-y_1)/uy;
        if (x == fmin(x,y)) {
            x0 = fmax(fx_1,fx_2);
            y0 = (x0-x_1)/ux*uy+y_1;
            z0 = (x0-x_1)/ux*uz+z_1;
        }
        else {
            y0 = fmax(fy_1, fy_2);
            x0 = (y0-y_1)/uy*ux+x_1;
            z0 = (y0-y_1)/uy*uz+z_1;
        }
    }
    else if ( (fy_1 != fy_2) &&(fz_1 != fz_2) ) { //#3
        fy_2=fy_1+SIGN(fy_2-fy_1);//added
        fz_2=fz_1+SIGN(fz_2-fz_1);//added
        y = (fmax(fy_1,fy_2)-y_1)/uy;
        z = (fmax(fz_1,fz_2)-z_1)/uz;
        if (y == fmin(y,z)) {
            y0 = fmax(fy_1,fy_2);
            x0 = (y0-y_1)/uy*ux+x_1;
            z0 = (y0-y_1)/uy*uz+z_1;
        }
        else {
            z0 = fmax(fz_1, fz_2);
            x0 = (z0-z_1)/uz*ux+x_1;
            y0 = (z0-z_1)/uz*uy+y_1;
        }
    }
    else if ( (fx_1 != fx_2) && (fz_1 != fz_2) ) { //#4
        fx_2=fx_1+SIGN(fx_2-fx_1);//added
        fz_2=fz_1+SIGN(fz_2-fz_1);//added
        x = (fmax(fx_1,fx_2)-x_1)/ux;
        z = (fmax(fz_1,fz_2)-z_1)/uz;
        if (x == fmin(x,z)) {
            x0 = fmax(fx_1,fx_2);
            y0 = (x0-x_1)/ux*uy+y_1;
            z0 = (x0-x_1)/ux*uz+z_1;
        }
        else {
            z0 = fmax(fz_1, fz_2);
            x0 = (z0-z_1)/uz*ux+x_1;
            y0 = (z0-z_1)/uz*uy+y_1;
        }
    }
    else if (fx_1 != fx_2) { //#5
        fx_2=fx_1+SIGN(fx_2-fx_1);//added
        x0 = fmax(fx_1,fx_2);
        y0 = (x0-x_1)/ux*uy+y_1;
        z0 = (x0-x_1)/ux*uz+z_1;
    }
    else if (fy_1 != fy_2) { //#6
        fy_2=fy_1+SIGN(fy_2-fy_1);//added
        y0 = fmax(fy_1, fy_2);
        x0 = (y0-y_1)/uy*ux+x_1;
        z0 = (y0-y_1)/uy*uz+z_1;
    }
    else { //#7 
        z0 = fmax(fz_1, fz_2);
        fz_2=fz_1+SIGN(fz_2-fz_1);//added
        x0 = (z0-z_1)/uz*ux+x_1;
        y0 = (z0-z_1)/uz*uy+y_1;
    }
    //s = ( (x0-fx_1)*dx + (y0-fy_1)*dy + (z0-fz_1)*dz )/3;
    //s = sqrt( SQR((x0-x_1)*dx) + SQR((y0-y_1)*dy) + SQR((z0-z_1)*dz) );
    //s = sqrt(SQR(x0-x_1)*SQR(dx) + SQR(y0-y_1)*SQR(dy) + SQR(z0-z_1)*SQR(dz));
    s = sqrt( SQR((x0-x_1)*dx) + SQR((y0-y_1)*dy) + SQR((z0-z_1)*dz));
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

