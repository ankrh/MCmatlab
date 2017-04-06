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
	cc -o mcxyz_unix mcxyz.c
This yields the executable mcxyz_unix (you can use any name).
	
2. Create .mci file
Run the MATLAB file 
	maketissue.m
which as an example yields 
	skinvessel_H.mci = monte carlo input file, a text file that gomcxyz will read,
	skinvessel_T.bin = binary file holding T(y,x,z), which identifies the tissue type within each voxel by an integer. Also produces a drawing of the tissue structure, and the rays of the illuminating beam.

3. Run the program on UNIX:
	mcxyz_unix skinvessel
This command will cause mcxyz_unix to read the skinvessel_H.mci and skinvessel_T.bin files, and run. The output is 
	skinvessel_F.bin = a binary file holding fluence rate F(y,x,z).
	skinvessel_props.m = a .m MATLAB file which loads the optical properties for the tissue types used in the example, muav(:),musv(:), gv(:). Just for convenience if one wishes to double-check the run parameters.

4. Look at results using MATLAB:
Run the MATLAB file
	lookmcxyz.m
which produces side views of the fluence rate, Fzx(z,x)@y=Ny/2 and Fzy(z,y)@x=Nx/2. Also produces a drawing of the tissue structure, and the rays of the illuminating beam.

	