# README #

### What is this repository for? ###

A Monte Carlo simulation for modeling light propagation in a 3D voxel space.
Fluorescence can optionally be simulated after simulation of the excitation light.
Included is also a finite element simulation for temperature increase and heat diffusion in the same voxel space.
Primarily targeted for tissue optics, but can be used in any environment.

### LICENSE ###

Attribution 4.0 International (CC BY 4.0)
You are welcome to use this code in whatever way you want. When doing so, please cite this code's doi: 10.13140/RG.2.2.32610.43205

### How do I get set up? ###

A quick description of how to run mcxyz, using MATLAB R2014B or later.
Note that if you want to use the pre-compiled mex-files, you may need to update MATLAB to the most recent version to avoid crashes.

PROGRAMS:
makeTissue.m
runMonteCarlo.m
runMonteCarloFluorescence.m
lookmcxyz.m
simulateHeatDistribution.m

HELPERFILE:
makeTissueList.m
	
INSTRUCTIONS:
1. Compilation
The folders include all the executables necessary, so you don't need to compile anything. If, however, you want to change the routine in either the mcxyz.c source code (in the folder "src") or the finiteElementHeatPropagator.c source code (in the folder "src"), you will need to recompile the respective mex-files. Check out those two source-files on how to do so.

2. Make Tissue
- Modify makeTissueList.m to include the definitions of the tissue types you're interested in. You may optionally include thermal and/or fluorescence properties.
- Modify makeTissue.m, to build up your tissue structure from the tissues defined in makeTissueList in the voxel space.
- Run makeTissue.m in MATLAB. This yields a .mat file named as defined in makeTissue.m (in the folder "Data"), containing the tissue definition.

3. Calculate light distribution
- Modify runMonteCarlo.m to include the beam type you're interested in, as well as the time you want to simulate photons. 1 minute should be sufficient for testing purposes. Find more instructions on how to define the beam in the file directly.
- Type "runMonteCarlo('[name]')" in the command prompt, with [name] being the name defined in makeTissue.m in step 2. This command will simulate photons for the requested time. The output is "name_MCoutput.mat" in the folder "Data", holding the relative fluence rate F(x,y,z).

4. Look at results using MATLAB:
- lookmcxyz('[name]') will be automatically called after runMonteCarlo, but can also be called independently to inspect the outcome of earlier Monte Carlo runs.

5. (Optional) Calculate fluorescence light distribution
- Modify runMonteCarloFluorescence.m to model fluorescence for a given input power.
- Type "runMonteCarloFluorescence('[name]')" in the command prompt, with [name] being the name defined in makeTissue.m in step 2. This command will simulate photons for the requested time. The output is "name_MCoutput_fluorescence.mat" in the folder "Data", holding the absolute fluorescence fluence rate I(x,y,z).

6. (Optional) Simulate Heat Distribution
Note: Heat simulations do not include fluorescence effects!
- Modify simulateHeatDistribution.m to specify the amount of illumination time (where energy is deposited as well as diffused) and dark time (where the lights are off, and the heat is just diffused). You can also specify how often you want a visual update of the simulation. You can also tell the script to generate a video (.mp4 format) of the temperature evolution. If the simulations are very slow, go back to makeTissue.m and try defining your voxel space with a smaller cuboid and/or lower resolution. 
- Type "[temperatureSensor, timeVector] = simulateHeatDistribution('[name]');" in the command prompt. A figure will appear, giving you the possibility to put one or more (shift-click) temperature sensors into the tissue, where the temperature throughout the simulation will be recorded and stored in the array temperatureSensor. After placing all required sensors, press any key. The simulation will then start running, updating the temperature distribution visualisation as often as you requested.

7. Look at all the results and start playing around.

### Contribution guidelines ###

This software is in productive use. When coding, make sure you don't brake the existing code (or, even better, fix the wreck), comment in the files what you changed, and commit.

### Who do I talk to? ###

This repository is part of the DTU Fotonik "biophotonics" team.
The main responsible is Anders Kragh Hansen: ankrh@fotonik.dtu.dk
