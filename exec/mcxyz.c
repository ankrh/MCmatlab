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
 **********/

#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>

#define Ntiss		19          /* Number of tissue types. */
#define STRLEN 		32          /* String length. */
#define ls          1.0E-7      /* Moving photon a little bit off the voxel face */
#define	PI          3.1415926
#define	LIGHTSPEED	2.997925E10 /* in vacuo speed of light [cm/s] */
#define ALIVE       1   		/* if photon not yet terminated */
#define DEAD        0    		/* if photon is to be terminated */
#define THRESHOLD   0.01		/* used in roulette */
#define CHANCE      0.1  		/* used in roulette */
#define Boolean     char
#define SQR(x)		(x*x) 
#define SIGN(x)     ((x)>=0 ? 1:-1)
#define RandomNum   (double) RandomGen(1, 0, NULL) /* Calls for a random number. */
#define COS90D      1.0E-6          /* If cos(theta) <= COS90D, theta >= PI/2 - 1e-6 rad. */
#define ONE_MINUS_COSZERO 1.0E-12   /* If 1-cos(theta) <= ONE_MINUS_COSZERO, fabs(theta) <= 1e-6 rad. */
/* If 1+cos(theta) <= ONE_MINUS_COSZERO, fabs(PI-theta) <= 1e-6 rad. */

/* DECLARE FUNCTIONS */
double RandomGen(char Type, long Seed, long *Status);  
/* Random number generator */
Boolean SameVoxel(double x1,double y1,double z1, double x2, double y2, double z2, double dx,double dy,double dz);
/* Asks,"In the same voxel?" */
double max2(double a, double b);
double min2(double a, double b);
double min3(double a, double b, double c);
double FindVoxelFace(double x1,double y1,double z1, double x2, double y2, double z2,double dx,double dy,double dz, double ux, double uy, double uz);
double FindVoxelFace2(double x1,double y1,double z1, double x2, double y2, double z2,double dx,double dy,double dz, double ux, double uy, double uz);
/* How much step size will the photon take to get the first voxel crossing in one single long step? */
double RFresnel(double n1, double n2, double ca1, double *ca2_Ptr);

int main(int argc, const char * argv[]) {
    
    if (argc==0) {
        printf("assuming you've compiled mcxyz.c as gomcxyz ...\n");
		printf("USAGE: gomcxyz name\n");
		printf("which will load the files name_H.mci and name_T.bin\n");
		printf("and run the Monte Carlo program.\n");
		printf("Yields  name_F.bin, which holds the fluence rate distribution.\n");
        return 0;
    }
	
	/* Propagation parameters */
	double	x, y, z;        /* photon position */
	double	ux, uy, uz;     /* photon trajectory as cosines */
	double  uxx, uyy, uzz;	/* temporary values used during SPIN */
	double	s;              /* step sizes. s = -log(RND)/mus [cm] */
	double  sleft;          /* dimensionless */
	double	costheta;       /* cos(theta) */
	double  sintheta;       /* sin(theta) */
	double	cospsi;         /* cos(psi) */
	double  sinpsi;         /* sin(psi) */
	double	psi;            /* azimuthal angle */
	long	i_photon;       /* current photon */
	double	W;              /* photon weight */
	double	absorb;         /* weighted deposited in a step due to absorption */
	short   photon_status;  /* flag = ALIVE=1 or DEAD=0 */
	Boolean sv;             /* Are they in the same voxel? */
	
	/* other variables */
	double	mua;            /* absorption coefficient [cm^-1] */
	double	mus;            /* scattering coefficient [cm^-1] */
	double	g;              /* anisotropy [-] */
	double	Nphotons;       /* number of photons in simulation */
	
	/* launch parameters */
	int		mcflag, launchflag, boundaryflag;
	float	xfocus, yfocus, zfocus;
	float	ux0, uy0, uz0;
	float	radius;
	float	waist;
	
	/* dummy variables */
	double  rnd;            /* assigned random value 0-1 */
	double	r, phi;			/* dummy values */
	long	i,j,NN;         /* dummy indices */
	double	tempx, tempy, tempz; /* temporary variables, used during photon step. */
	int 	ix, iy, iz;     /* Added. Used to track photons */
	double 	temp;           /* dummy variable */
    int     bflag;          /* boundary flag:  0 = photon inside volume. 1 = outside volume */
	int		CNT;
	
	/* mcxyz bin variables */
	float	dx, dy, dz;     /* bin size [cm] */
	int		Nx, Ny, Nz, Nt; /* # of bins */
	float	xs, ys, zs;		/* launch position */
    
    /* time */
	float	time_min;               // Requested time duration of computation.
	time_t	now;
	double	start_time, finish_time, temp_time; /* for clock() */
	
	/* tissue parameters */
	char	tissuename[50][32];
	float 	muav[Ntiss];            // muav[0:Ntiss-1], absorption coefficient of ith tissue type
	float 	musv[Ntiss];            // scattering coeff. 
	float 	gv[Ntiss];              // anisotropy of scattering
    
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
		sscanf(buf, "%f", &time_min); // desired time duration of run [min]
		fgets(buf, 32, fid);
		sscanf(buf, "%d", &Nx);  // # of bins  
		fgets(buf, 32,fid);
		sscanf(buf, "%d", &Ny);  // # of bins
		fgets(buf, 32,fid);
		sscanf(buf, "%d", &Nz);  // # of bins   
	
		fgets(buf, 32,fid);
		sscanf(buf, "%f", &dx);	 // size of bins [cm]
		fgets(buf, 32,fid);
		sscanf(buf, "%f", &dy);	 // size of bins [cm] 
		fgets(buf, 32,fid);
		sscanf(buf, "%f", &dz);	 // size of bins [cm] 
	
		// launch parameters
		fgets(buf, 32,fid);
		sscanf(buf, "%d", &mcflag);  // mcflag, 0 = uniform, 1 = Gaussian, 2 = iso-pt
		fgets(buf, 32,fid);
		sscanf(buf, "%d", &launchflag);  // launchflag, 0 = ignore, 1 = manually set
        fgets(buf, 32,fid);
        sscanf(buf, "%d", &boundaryflag);  // 0 = no boundaries, 1 = escape at all boundaries, 2 = escape at surface only
	
		fgets(buf, 32,fid);
		sscanf(buf, "%f", &xs);  // initial launch point
		fgets(buf, 32,fid);
		sscanf(buf, "%f", &ys);  // initial launch point 
		fgets(buf, 32,fid);
		sscanf(buf, "%f", &zs);  // initial launch point
	
		fgets(buf, 32,fid);
		sscanf(buf, "%f", &xfocus);  // xfocus
		fgets(buf, 32,fid);
		sscanf(buf, "%f", &yfocus);  // yfocus
		fgets(buf, 32,fid);
		sscanf(buf, "%f", &zfocus);  // zfocus

		fgets(buf, 32,fid);
		sscanf(buf, "%f", &ux0);  // ux trajectory
		fgets(buf, 32,fid);
		sscanf(buf, "%f", &uy0);  // uy trajectory
		fgets(buf, 32,fid);
		sscanf(buf, "%f", &uz0);  // uz trajectory

		fgets(buf, 32,fid);
		sscanf(buf, "%f", &radius);  // radius
		fgets(buf, 32,fid);
		sscanf(buf, "%f", &waist);  // waist
    
	
		// tissue optical properties
		fgets(buf, 32,fid);
		sscanf(buf, "%d", &Nt);				// # of tissue types in tissue list
		for (i=1; i<=Nt; i++) {
			fgets(buf, 32, fid);
			sscanf(buf, "%f", &muav[i]);	// absorption coeff [cm^-1]
			fgets(buf, 32, fid);
			sscanf(buf, "%f", &musv[i]);	// scattering coeff [cm^-1]
			fgets(buf, 32, fid);
			sscanf(buf, "%f", &gv[i]);		// anisotropy of scatter [dimensionless]
		}    
    fclose(fid);
    
    printf("time_min = %0.2f min\n",time_min);
    printf("Nx = %d, dx = %0.4f [cm]\n",Nx,dx);
    printf("Ny = %d, dy = %0.4f [cm]\n",Ny,dy);
    printf("Nz = %d, dz = %0.4f [cm]\n",Nz,dz);

    printf("xs = %0.4f [cm]\n",xs);
    printf("ys = %0.4f [cm]\n",ys);
    printf("zs = %0.4f [cm]\n",zs);
    printf("mcflag = %d [cm]\n",mcflag);
    if (mcflag==0) printf("launching uniform flat-field beam\n");
    if (mcflag==1) printf("launching Gaissian beam\n");
    if (mcflag==2) printf("launching isotropic point source\n");
    printf("xfocus = %0.4f [cm]\n",xfocus);
    printf("yfocus = %0.4f [cm]\n",yfocus);
    printf("zfocus = %0.2e [cm]\n",zfocus);
	if (launchflag==1) {
		printf("Launchflag ON, so launch the following:\n");
		printf("ux0 = %0.4f [cm]\n",ux0);
		printf("uy0 = %0.4f [cm]\n",uy0);
		printf("uz0 = %0.4f [cm]\n",uz0);
	}
	else {
		printf("Launchflag OFF, so program calculates launch angles.\n");
		printf("radius = %0.4f [cm]\n",radius);
		printf("waist  = %0.4f [cm]\n",waist);
	}
    if (boundaryflag==0)
		printf("boundaryflag = 0, so no boundaries.\n");
    else if (boundaryflag==1)
		printf("boundaryflag = 1, so escape at all boundaries.\n");    
	else if (boundaryflag==2)
		printf("boundaryflag = 2, so escape at surface only.\n");    
	else{
        printf("improper boundaryflag. quit.\n");
        return 0;
    }
    printf("# of tissues available, Nt = %d\n",Nt);
    for (i=1; i<=Nt; i++) {
        printf("muav[%ld] = %0.4f [cm^-1]\n",i,muav[i]);
        printf("musv[%ld] = %0.4f [cm^-1]\n",i,musv[i]);
        printf("  gv[%ld] = %0.4f [--]\n\n",i,gv[i]);
    }

    // SAVE optical properties, for later use by MATLAB.
	strcpy(filename,myname);
	strcat(filename,"_props.m");
	fid = fopen(filename,"w");
	for (i=1; i<=Nt; i++) {
		fprintf(fid,"muav(%ld) = %0.4f;\n",i,muav[i]);
		fprintf(fid,"musv(%ld) = %0.4f;\n",i,musv[i]);
		fprintf(fid,"gv(%ld) = %0.4f;\n\n",i,gv[i]);
	}
	fclose(fid);
    
    /* IMPORT BINARY TISSUE FILE */
	char 	*v=NULL;
	float 	*F=NULL;
	float	*R=NULL;
	int 	type;
	NN = Nx*Ny*Nz;
	v  = ( char *)malloc(NN*sizeof(char));  /* tissue structure */
	F  = (float *)malloc(NN*sizeof(float));	/* relative fluence rate [W/cm^2/W.delivered] */
	// NOT READY: R  = (float *)malloc(NN*sizeof(float));	/* escaping flux [W/cm^2/W.delivered] */
    
	// read binary file
    strcpy(filename,myname);
    strcat(filename, "_T.bin");
    fid = fopen(filename, "rb");
    fread(v, sizeof(char), NN, fid);
    fclose(fid);
    
    // Show tissue on screen, along central z-axis, by listing tissue type #'s.
    iy = Ny/2;
    ix = Nx/2;
    printf("central axial profile of tissue types:\n");
    for (iz=0; iz<Nz; iz++) {
        i = (long)(iz*Ny*Nx + ix*Ny + iy);
        printf("%d",v[i]);
    }
    printf("\n\n");
    
	/**************************
	 * ============================ MAJOR CYCLE ========================
	 **********/
	start_time = clock();
	now = time(NULL);
	printf("%s\n", ctime(&now));	
	
	/**** INITIALIZATIONS 
	 *****/
	RandomGen(0, -(int)time(NULL)%(1<<15), NULL); /* initiate with seed = 1, or any long integer. */
	for(j=0; j<NN;j++) 	F[j] = 0; // ensure F[] starts empty.	
	
	/**** RUN
	 Launch N photons, initializing each one before progation.
	 *****/
	printf("------------- Begin Monte Carlo -------------\n");
    printf("%s\n",myname);
    printf("requesting %0.1f min\n",time_min);
	Nphotons = 1000; // will be updated to achieve desired run time, time_min.
	i_photon = 0;
	CNT = 0;
	do {
		/**** LAUNCH 
		 Initialize photon position and trajectory.
		 *****/
		//if (fmod(i_photon,10)==0) printf("photon %ld took %d steps\n",i_photon,CNT);

		i_photon += 1;				/* increment photon count */
		W = 1.0;                    /* set photon weight to one */
		photon_status = ALIVE;      /* Launch an ALIVE photon */
		CNT = 0;
		
		// Print out message about progress.
		if ((i_photon>1000) & (fmod(i_photon, (int)(Nphotons/100))  == 0)) {
            temp = i_photon/Nphotons*100;
            //printf("%0.1f%% \t\tfmod = %0.3f\n", temp,fmod(temp, 10.0));
            if ((temp<10) | (temp>90)){
                printf("%0.0f%% done\n", i_photon/Nphotons*100);
            }
            else if(fmod(temp, 10.0)>9)
                printf("%0.0f%% done\n", i_photon/Nphotons*100);
        }
        
		// At 1000th photon, update Nphotons to achieve desired runtime (time_min)
		if (i_photon==1)
			temp_time = clock();
		if (i_photon==1000) {    
			finish_time = clock();
			Nphotons = (long)( time_min*60*999*CLOCKS_PER_SEC/(finish_time-temp_time) );
			printf("Nphotons = %0.0f for simulation time = %0.2f min\n",Nphotons,time_min);
		}
		
		/**** SET SOURCE 
		 * Launch collimated beam at x,y center. 
		 ****/
        
		/****************************/
		/* Initial position. */		
		
		/* trajectory */
		if (launchflag==1) { // manually set launch
			x	= xs; 
			y	= ys;
			z	= zs;
			ux	= ux0;
			uy	= uy0;
			uz	= uz0;
		}
		else { // use mcflag
			if (mcflag==0) { // uniform beam
				// set launch point and width of beam
				while ((rnd = RandomGen(1,0,NULL)) <= 0.0); // avoids rnd = 0
				r		= radius*sqrt(rnd); // radius of beam at launch point
				while ((rnd = RandomGen(1,0,NULL)) <= 0.0); // avoids rnd = 0
				phi		= rnd*2.0*PI;
				x		= xs + r*cos(phi);
				y		= ys + r*sin(phi);
				z		= zs;
				// set trajectory toward focus
				while ((rnd = RandomGen(1,0,NULL)) <= 0.0); // avoids rnd = 0
				r		= waist*sqrt(rnd); // radius of beam at focus
				while ((rnd = RandomGen(1,0,NULL)) <= 0.0); // avoids rnd = 0
				phi		= rnd*2.0*PI;
				xfocus	= r*cos(phi);
				yfocus	= r*sin(phi);
				temp	= sqrt((x - xfocus)*(x - xfocus) + (y - yfocus)*(y - yfocus) + zfocus*zfocus);
				ux		= -(x - xfocus)/temp;
				uy		= -(y - yfocus)/temp;
				uz		= sqrt(1 - ux*ux + uy*uy);
			}
			else { // isotropic pt source
				costheta = 1.0 - 2.0*RandomGen(1,0,NULL);
				sintheta = sqrt(1.0 - costheta*costheta);
				psi = 2.0*PI*RandomGen(1,0,NULL);
				cospsi = cos(psi);
				if (psi < PI)
					sinpsi = sqrt(1.0 - cospsi*cospsi); 
				else
					sinpsi = -sqrt(1.0 - cospsi*cospsi);
				x = xs;
				y = ys;
				z = zs;
				ux = sintheta*cospsi;
				uy = sintheta*sinpsi;
				uz = costheta;
			}
		} // end  use mcflag
		/****************************/
		
		/* Get tissue voxel properties of launchpoint.
		 * If photon beyond outer edge of defined voxels, 
		 * the tissue equals properties of outermost voxels.
		 * Therefore, set outermost voxels to infinite background value.
		 */
		ix = (int)(Nx/2 + x/dx);
		iy = (int)(Ny/2 + y/dy);
		iz = (int)(z/dz);        
		if (ix>=Nx) ix=Nx-1;
		if (iy>=Ny) iy=Ny-1;
		if (iz>=Nz) iz=Nz-1;
		if (ix<0)   ix=0;
		if (iy<0)   iy=0;
		if (iz<0)   iz=0;		
		/* Get the tissue type of located voxel */
		i		= (long)(iz*Ny*Nx + ix*Ny + iy);
		type	= v[i];
		mua 	= muav[type];
		mus 	= musv[type];
		g 		= gv[type];
		
        bflag = 1; // initialize as 1 = inside volume, but later check as photon propagates.
        
		/* HOP_DROP_SPIN_CHECK
		 Propagate one photon until it dies as determined by ROULETTE.
		 *******/
		do {
			
			/**** HOP
			 Take step to new position
			 s = dimensionless stepsize
			 x, uy, uz are cosines of current photon trajectory
			 *****/
			while ((rnd = RandomNum) <= 0.0);   /* yields 0 < rnd <= 1 */
			sleft	= -log(rnd);				/* dimensionless step */
			CNT += 1;
			
			do{  // while sleft>0   
				s     = sleft/mus;				/* Step size [cm].*/
				tempx = x + s*ux;				/* Update positions. [cm] */
				tempy = y + s*uy;	
				tempz = z + s*uz;
				
				sv = SameVoxel(x,y,z, tempx, tempy, tempz, dx,dy,dz);
				if (sv) /* photon in same voxel */
				{  
					x=tempx;					/* Update positions. */
					y=tempy;
					z=tempz;
					
					/**** DROP
					 Drop photon weight (W) into local bin.
					 *****/
                    absorb = W*(1 - exp(-mua*s));	/* photon weight absorbed at this step */
                    W -= absorb;					/* decrement WEIGHT by amount absorbed */
					// If photon within volume of heterogeneity, deposit energy in F[]. 
					// Normalize F[] later, when save output. 
                    if (bflag) F[i] += absorb;	// only save data if blag==1, i.e., photon inside simulation cube
					
					/* Update sleft */
					sleft = 0;		/* dimensionless step remaining */
				}
				else /* photon has crossed voxel boundary */
				{
					/* step to voxel face + "littlest step" so just inside new voxel. */
					s = ls + FindVoxelFace2(x,y,z, tempx,tempy,tempz, dx,dy,dz, ux,uy,uz);
					
					/**** DROP
					 Drop photon weight (W) into local bin.
					 *****/
					absorb = W*(1-exp(-mua*s));   /* photon weight absorbed at this step */
					W -= absorb;                  /* decrement WEIGHT by amount absorbed */
					// If photon within volume of heterogeneity, deposit energy in F[]. 
					// Normalize F[] later, when save output. 
                    if (bflag) F[i] += absorb;	
					
					/* Update sleft */
					sleft -= s*mus;  /* dimensionless step remaining */
					if (sleft<=ls) sleft = 0;
					
					/* Update positions. */
					x += s*ux;
					y += s*uy;
					z += s*uz;
					
					// pointers to voxel containing optical properties
                    ix = (int)(Nx/2 + x/dx);
                    iy = (int)(Ny/2 + y/dy);
                    iz = (int)(z/dz);
                    
                    bflag = 1;  // Boundary flag. Initialize as 1 = inside volume, then check.
                    if (boundaryflag==0) { // Infinite medium.
								// Check if photon has wandered outside volume.
                                // If so, set tissue type to boundary value, but let photon wander.
                                // Set blag to zero, so DROP does not deposit energy.
						if (iz>=Nz) {iz=Nz-1; bflag = 0;}
						if (ix>=Nx) {ix=Nx-1; bflag = 0;}
						if (iy>=Ny) {iy=Ny-1; bflag = 0;}
						if (iz<0)   {iz=0;    bflag = 0;}
						if (ix<0)   {ix=0;    bflag = 0;}
						if (iy<0)   {iy=0;    bflag = 0;}
					}
					else if (boundaryflag==1) { // Escape at boundaries
						if (iz>=Nz) {iz=Nz-1; photon_status = DEAD; sleft=0;}
						if (ix>=Nx) {ix=Nx-1; photon_status = DEAD; sleft=0;}
						if (iy>=Ny) {iy=Ny-1; photon_status = DEAD; sleft=0;}
						if (iz<0)   {iz=0;    photon_status = DEAD; sleft=0;}
						if (ix<0)   {ix=0;    photon_status = DEAD; sleft=0;}
						if (iy<0)   {iy=0;    photon_status = DEAD; sleft=0;}
					}
					else if (boundaryflag==2) { // Escape at top surface, no x,y bottom z boundaries
						if (iz>=Nz) {iz=Nz-1; bflag = 0;}
						if (ix>=Nx) {ix=Nx-1; bflag = 0;}
						if (iy>=Ny) {iy=Ny-1; bflag = 0;}
						if (iz<0)   {iz=0;    photon_status = DEAD; sleft=0;}
						if (ix<0)   {ix=0;    bflag = 0;}
						if (iy<0)   {iy=0;    bflag = 0;}
					}
					
                    // update pointer to tissue type
					i    = (long)(iz*Ny*Nx + ix*Ny + iy);
					type = v[i];
                    mua  = muav[type];
                    mus  = musv[type];
                    g    = gv[type];
                    
				} //(sv) /* same voxel */
                
			} while(sleft>0); //do...while
			
			/**** SPIN 
			 Scatter photon into new trajectory defined by theta and psi.
			 Theta is specified by cos(theta), which is determined 
			 based on the Henyey-Greenstein scattering function.
			 Convert theta and psi into cosines ux, uy, uz. 
			 *****/
			/* Sample for costheta */
			rnd = RandomNum;
			if (g == 0.0)
				costheta = 2.0*rnd - 1.0;
			else {
				double temp = (1.0 - g*g)/(1.0 - g + 2*g*rnd);
				costheta = (1.0 + g*g - temp*temp)/(2.0*g);
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
				temp = sqrt(1.0 - uz * uz);
				uxx = sintheta * (ux * uz * cospsi - uy * sinpsi) / temp + ux * costheta;
				uyy = sintheta * (uy * uz * cospsi + ux * sinpsi) / temp + uy * costheta;
				uzz = -sintheta * cospsi * temp + uz * costheta;
			}
			
			/* Update trajectory */
			ux = uxx;
			uy = uyy;
			uz = uzz;
			
			/**** CHECK ROULETTE 
			 If photon weight below THRESHOLD, then terminate photon using Roulette technique.
			 Photon has CHANCE probability of having its weight increased by factor of 1/CHANCE,
			 and 1-CHANCE probability of terminating.
			 *****/
			if (W < THRESHOLD) {
				if (RandomNum <= CHANCE)
					W /= CHANCE;
				else photon_status = DEAD;
			}				
            
		} while (photon_status == ALIVE);  /* end STEP_CHECK_HOP_SPIN */
        /* if ALIVE, continue propagating */
		/* If photon DEAD, then launch new photon. */	
        
	} while (i_photon < Nphotons);  /* end RUN */
	
    
	printf("------------------------------------------------------\n");
	finish_time = clock();
	time_min = (double)(finish_time-start_time)/CLOCKS_PER_SEC/60;
	printf("Elapsed Time for %0.3e photons = %5.3f min\n",Nphotons,time_min);
	printf("%0.2e photons per minute\n", Nphotons/time_min);
	
    /**** SAVE
     Convert data to relative fluence rate [cm^-2] and save.
     *****/
    
    // Normalize deposition (A) to yield fluence rate (F).
    temp = dx*dy*dz*Nphotons;
    for (i=0; i<NN;i++){
        F[i] /= (temp*muav[v[i]]);
    }
    // Save the binary file
    strcpy(filename,myname);
    strcat(filename,"_F.bin");
    printf("saving %s\n",filename);
    fid = fopen(filename, "wb");   /* 3D voxel output */
    fwrite(F, sizeof(float), NN, fid);
    fclose(fid);
    
    /* save reflectance */
// NOT READY: 
//	strcpy(filename,myname);
//    strcat(filename,"_Ryx.bin");
//    printf("saving %s\n",filename);
//    fid = fopen(filename, "wb");   /* 2D voxel output */
//	int Nyx = Ny*Nx;
//    fwrite(R, sizeof(float), Nyx, fid);
//    fclose(fid);
//    printf("%s is done.\n",myname);
	
    printf("------------------------------------------------------\n");
    now = time(NULL);
    printf("%s\n", ctime(&now));
    
    
    free(v);
    free(F);
    free(R);
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
#define MBIG 1000000000
#define MSEED 161803398
#define MZ 0
#define FAC 1.0E-9

double RandomGen(char Type, long Seed, long *Status){
    static long i1, i2, ma[56];   /* ma[0] is not used. */
    long        mj, mk;
    short       i, ii;
    
    if (Type == 0) {              /* set seed. */
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
Boolean SameVoxel(double x1,double y1,double z1, double x2, double y2, double z2, double dx,double dy,double dz)
{
    double xmin=min2((floor)(x1/dx),(floor)(x2/dx))*dx;
    double ymin=min2((floor)(y1/dy),(floor)(y2/dy))*dy;
    double zmin=min2((floor)(z1/dz),(floor)(z2/dz))*dz;
    double xmax = xmin+dx;
    double ymax = ymin+dy;
    double zmax = zmin+dz;
    Boolean sv=0;
    
    sv=(x1<=xmax && x2<=xmax && y1<=ymax && y2<=ymax && z1<zmax && z2<=zmax);
    return (sv);
}

/***********************************************************
 * max2
 ****/
double max2(double a, double b) {
    double m;
    if (a > b)
        m = a;
    else
        m = b;
    return m;
}

/***********************************************************
 * min2
 ****/
double min2(double a, double b) {
    double m;
    if (a >= b)
        m = b;
    else
        m = a;
    return m;
}
/***********************************************************
 * min3
 ****/
double min3(double a, double b, double c) {
    double m;
    if (a <=  min2(b, c))
        m = a;
    else if (b <= min2(a, c))
        m = b;
    else
        m = c;
    return m;
}

/********************
 * my version of FindVoxelFace for no scattering.
 * s = ls + FindVoxelFace2(x,y,z, tempx, tempy, tempz, dx, dy, dz, ux, uy, uz);
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
        
        x = (max2(fx_1,fx_2)-x_1)/ux;
        y = (max2(fy_1,fy_2)-y_1)/uy;
        z = (max2(fz_1,fz_2)-z_1)/uz;
        if (x == min3(x,y,z)) {
            x0 = max2(fx_1,fx_2);
            y0 = (x0-x_1)/ux*uy+y_1;
            z0 = (x0-x_1)/ux*uz+z_1;
        }
        else if (y == min3(x,y,z)) {
            y0 = max2(fy_1,fy_2);
            x0 = (y0-y_1)/uy*ux+x_1;
            z0 = (y0-y_1)/uy*uz+z_1;
        }
        else {
            z0 = max2(fz_1,fz_2);
            y0 = (z0-z_1)/uz*uy+y_1;
            x0 = (z0-z_1)/uz*ux+x_1;
        }
    }
    else if ( (fx_1 != fx_2) && (fy_1 != fy_2) ) { //#2
        fx_2=fx_1+SIGN(fx_2-fx_1);//added
        fy_2=fy_1+SIGN(fy_2-fy_1);//added
        x = (max2(fx_1,fx_2)-x_1)/ux;
        y = (max2(fy_1,fy_2)-y_1)/uy;
        if (x == min2(x,y)) {
            x0 = max2(fx_1,fx_2);
            y0 = (x0-x_1)/ux*uy+y_1;
            z0 = (x0-x_1)/ux*uz+z_1;
        }
        else {
            y0 = max2(fy_1, fy_2);
            x0 = (y0-y_1)/uy*ux+x_1;
            z0 = (y0-y_1)/uy*uz+z_1;
        }
    }
    else if ( (fy_1 != fy_2) &&(fz_1 != fz_2) ) { //#3
        fy_2=fy_1+SIGN(fy_2-fy_1);//added
        fz_2=fz_1+SIGN(fz_2-fz_1);//added
        y = (max2(fy_1,fy_2)-y_1)/uy;
        z = (max2(fz_1,fz_2)-z_1)/uz;
        if (y == min2(y,z)) {
            y0 = max2(fy_1,fy_2);
            x0 = (y0-y_1)/uy*ux+x_1;
            z0 = (y0-y_1)/uy*uz+z_1;
        }
        else {
            z0 = max2(fz_1, fz_2);
            x0 = (z0-z_1)/uz*ux+x_1;
            y0 = (z0-z_1)/uz*uy+y_1;
        }
    }
    else if ( (fx_1 != fx_2) && (fz_1 != fz_2) ) { //#4
        fx_2=fx_1+SIGN(fx_2-fx_1);//added
        fz_2=fz_1+SIGN(fz_2-fz_1);//added
        x = (max2(fx_1,fx_2)-x_1)/ux;
        z = (max2(fz_1,fz_2)-z_1)/uz;
        if (x == min2(x,z)) {
            x0 = max2(fx_1,fx_2);
            y0 = (x0-x_1)/ux*uy+y_1;
            z0 = (x0-x_1)/ux*uz+z_1;
        }
        else {
            z0 = max2(fz_1, fz_2);
            x0 = (z0-z_1)/uz*ux+x_1;
            y0 = (z0-z_1)/uz*uy+y_1;
        }
    }
    else if (fx_1 != fx_2) { //#5
        fx_2=fx_1+SIGN(fx_2-fx_1);//added
        x0 = max2(fx_1,fx_2);
        y0 = (x0-x_1)/ux*uy+y_1;
        z0 = (x0-x_1)/ux*uz+z_1;
    }
    else if (fy_1 != fy_2) { //#6
        fy_2=fy_1+SIGN(fy_2-fy_1);//added
        y0 = max2(fy_1, fy_2);
        x0 = (y0-y_1)/uy*ux+x_1;
        z0 = (y0-y_1)/uy*uz+z_1;
    }
    else { //#7 
        z0 = max2(fz_1, fz_2);
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