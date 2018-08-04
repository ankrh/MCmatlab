# README #

### What is this repository for? ###

A Monte Carlo simulation for modeling light propagation in a 3D voxel space.
Fluorescence can optionally be simulated after simulation of the excitation light.
Included is also a finite element simulation for temperature increase and heat diffusion in the same voxel space.
Primarily targeted for tissue optics, but can be used in any environment.

### LICENSE ###

TBD

### How do I get set up? ###

Requirements:
- Windows 7 or later, or macOS 10.12 (Sierra) or later (this software is not compatible with macOS 10.11 El Capitan)
- MATLAB R2014B or later (if you want to use the pre-compiled mex-files, you need to update MATLAB to the most recent version)
 - Image Processing Toolbox for MATLAB

MATLAB PROGRAMS:
defineGeometry.m
runMonteCarlo.m
runMonteCarloFluorescence.m
plotMCmatlab.m
simulateHeatDistribution.m

HELPERFILE:
getMediaProperties.m

INSTRUCTIONS:
1. Compilation
The folders include all the executables necessary, so you don't need to compile anything. If, however, you want to change the routine in either the MCmatlab.c source code or the finiteElementHeatPropagator.c source code (both located in the folder "src"), you will need to recompile the respective mex-files. Check out those two source-files on how to do so.

2. Build geometry
 - Modify getMediaProperties.m to include the definitions of the media you're interested in. You may optionally include thermal and/or fluorescence properties.
 - Modify defineGeometry.m, to build up your geometry from the media defined in getMediaProperties in the voxel space.
 - Run defineGeometry.m in MATLAB. This yields a .mat file named as defined in defineGeometry.m (in the folder "Data"), containing the geometry definition.

3. Calculate light distribution
 - Modify runMonteCarlo.m to include the beam type you're interested in, as well as the time you want to simulate photons. 0.1 minute should be sufficient for testing purposes. Find more instructions on how to define the beam in the file directly.
 - Type "runMonteCarlo('[name]')" in the command prompt, with [name] being the name defined in defineGeometry.m in step 2. This command will simulate photons for the requested time. The output is "name_MCoutput.mat" in the folder "Data", holding a struct containing the relative fluence rate F(x,y,z).

4. Look at results using MATLAB:
 - plotMCmatlab('[name]') will be automatically called after runMonteCarlo, but can also be called independently to inspect the outcome of earlier Monte Carlo runs.

5. (Optional) Calculate fluorescence light distribution
 - Modify runMonteCarloFluorescence.m to model fluorescence for a given input power.
 - Type "runMonteCarloFluorescence('[name]')" in the command prompt, with [name] being the name defined in defineGeometry.m in step 2. This command will simulate photons for the requested time. The output is "name_MCoutput_fluorescence.mat" in the folder "Data", holding a struct containing the absolute fluorescence fluence rate I(x,y,z).

6. (Optional, does not include fluorescence) Simulate heat distribution
 - Modify simulateHeatDistribution.m to specify the illumination time (where energy is deposited as well as diffused) and dark time (where the lights are off, and the heat is just diffused). You can also specify how often you want a visual update of the simulation. You can also tell the script to generate a video (.mp4 format) of the temperature evolution. Also, you can specify (x,y,z) locations for virtual temperature sensors in order to record, show and store the temperature as a function of time at those positions. If the simulations are very slow, go back to defineGeometry.m and try defining your voxel space with a smaller cuboid and/or lower resolution. 
 - Type "simulateHeatDistribution('[name]')" in the command prompt. The simulation will run, updating the temperature distribution visualisation as often as you requested.

7. Look at all the results and start playing around.

### Contribution guidelines ###

This software is in productive use. When coding, make sure you don't break the existing code (or, even better, fix the wreck), comment in the files what you changed, and commit.

### Who do I talk to? ###

This repository is part of the DTU Fotonik "biophotonics" team.
The main responsible is Anders Kragh Hansen: ankrh@fotonik.dtu.dk
