Readme.txt
A quick description of how to run mcxyz.c, 
using UNIX and MATLAB.

PROGRAMS:
makeTissueList.m
makeTissue.m
lookmcxyz.m
	
mcxyz.c

INSTRUCTIONS:
1. Compile mcxyz.c for UNIX
	cc -o mcxyz -O3 mcxyz.c
This yields the executable mcxyz (you can use any name).
	
2. Create .mci file
- Modify makeTissueList.m, to include the definitions of the tissue types you're interested in.
- Modify makeTissue.m, to build up your tissue structure in the voxel space. makeTissue.m also includes the definition of the beam you want to launch into that tissue voxel space.
- Run makeTissue.m in MATLAB. This, as an example, yields the following two files:
	skinvessel_H.mci = monte carlo input file, a text file that mcxyz will read,
	skinvessel_T.bin = binary file holding T(y,x,z), which identifies the tissue type within each voxel by an integer.

3. Make sure the two files generated in step 2 are in the same directory as mcxyz, and run the main program on a command line with the base name of your files as the first and only input:
	mcxyz skinvessel
This command will cause mcxyz to read the skinvessel_H.mci and skinvessel_T.bin files, and run. The output is 
	skinvessel_F.bin = a binary file holding fluence rate F(y,x,z).
	skinvessel_props.m = a .m MATLAB file which loads the optical properties for the tissue types used in the example, muav(:),musv(:), gv(:). Just for convenience if one wishes to double-check the run parameters.

4. Look at results using MATLAB:
Run the MATLAB file
	lookmcxyz.m
which produces volumetric views of the fluence rate, Fzx(z,x)@y=Ny/2 and Fzy(z,y)@x=Nx/2. This also produces a drawing of the tissue structure.

	