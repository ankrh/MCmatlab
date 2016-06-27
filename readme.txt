Readme.txt
A quick description of how to run mcxyz.c, 
using UNIX and MATLAB.

PROGRAMS:
maketissue.m
	SpectralLIB.mat
	makeTissueList.m
	makecmap.m
mcxyz.c
lookmcxyz.m
	reportHmci.m
	makec2f.m

INSTRUCTIONS:
1. Compile mcxyz.c for UNIX
	cc -o gomcxyz mcxyz.c
	This yields the executable gomcxyz (you can use any name).
	
2. Create .mci file
	Run the MATLAB file 
		maketissue.m
	which as an example yields 
		skinvessel_H.mci = monte carlo input file, a text file that gomcxyz will read,
		skinvessel_T.bin = binary file holding T(y,x,z), which identifies the tissue type
							within each voxel by an integer.
	Also produces a drawing of the tissue structure, and the rays of the illuminating beam.

3. Run the program on UNIX:
	gomcxyz kinvessel
	This command will cause gomcxyz to read the kinvessel_H.mci and kinvessel_T.bin files, and run.
	The output is 
		skinvessel_F.bin = binarly file holding fluence rate F(y,x,z).
		skinvessel_props.m = .m MATLAB file which loads the optical properties  
							for the tissue types used in the example, muav(:),musv(:), gv(:).

4. Look at results using MATLAB:
	Run the MATLAB file
		lookmcxyz.m
	which produces side views of the fluence rate, Fzx(z,x)@y=Ny/2 and Fzy(z,y)@x=Nx/2.
	Also produces a drawing of the tissue structure, and the rays of the illuminating beam.
		
	