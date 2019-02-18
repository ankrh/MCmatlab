# MCmatlab #

## What is this repository for? ##

A Monte Carlo simulation for modeling light propagation in a 3D voxel space.
Fluorescence can optionally be simulated after simulation of the excitation light.
Included is also a finite element simulation for temperature increase and heat diffusion in the same voxel space.
Primarily targeted for tissue optics, but can be used in any environment.

You can find the latest version of MCmatlab and more information in the wiki on https://gitlab.gbar.dtu.dk/biophotonics/MCmatlab.
If you publish results obtained with this software, we would be thankful if you cited its accompanying article: doi:10.1117/1.JBO.23.12.121622

## LICENSE ##

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.

## How do I get set up? ##

Requirements:
- Windows 7 or macOS 10.12 (Sierra) or newer
- MATLAB R2017a or newer

HELPER FILES:
All the helper functions needed for running MCmatlab are located in the folder "helperfuncs", which therefore has to be on your MATLAB path. You will not need to modify any of these files (except for "getMediaProperties.m" as mentioned below). The example model files automatically add "helperfuncs" to the MATLAB path as the first command, and you should keep that practice also in your own model files.

MODEL FILES:
In MCmatlab, you set up your model in a single m-file. You can find a few examples to get you started in the root folder of MCmatlab. Once you're familiar with those, you can check the model file "Template.m" for all the possible switches you might use, but it is itself not a valid model file, as it also contains some mutually exclusive switches.

MEDIA PROPERTIES:
The optical properties of all media are defined in the file "getMediaProperties.m" in the folder "helperfuncs". Many media are already defined therein, and you can either modify the properties of those, or add more media types to this file. Make sure that each media has its distinct "j"-number, as this is the number you will refer to when building your model.

## How do I use MCmatlab? ##
### Compilation ###
- The folders include all the executables necessary, so you don't need to compile anything. If, however, you want to change the routine in either the MCmatlab.c source code or the finiteElementHeatPropagator.c source code (both located in the folder "src"), you will need to recompile the respective mex-files. Check out those two source-files on how to do so.
 
### Building the model ###
You build a model in a seperate m-file. Each model requires the first two and optionally more of the following steps:
1. Build geometry (check out "Example1_StandardTissue.m", section "%% Geometry definition")
 - Modify or add to "helperfuncs/getMediaProperties.m" to include the definitions of the media you're interested in. You may optionally include thermal and/or fluorescence properties.
 - In your model file, specify the geometry of your model.
 - The "GeomFunc" you define in the model file simply defines a 3D-matrix containing the media definition for each voxel in the model cube. The media are referred to by a number, corresponding to the "j" defined in "helperfuncs/getMediaProperties.m".
 - After this step, you will be shown two figures with the geometry you defined and an overview of the optical properties.

2. Calculate light distribution (check out "Example1_StandardTissue.m", section "%% Monte Carlo simulation")
 - This section in the model file contains the definitions for the Monte Carlo simulation.
 - After this step, you will be shown two figures with the normalized fluence rate and the absorbed light in the cube.
 
3. (Optional) Include Fresnel reflection and refraction
 - Check out "Example2_RefractionReflection.m" on how to implement changing refractive indices in your model.
 
4. (Optional) Simulate heat distribution (check out "Example3_BloodVessel.m", section "%% Heat simulation")
 - To simulate the heat diffusion, all media involved in the simulation must have the thermal properties volumetric heat capacity (VHC) and thermal conductivity (TC) defined in "helperfuncs/getMediaProperties.m".
 - If you want to model chemical changes such as tissue damage, the media that might undergo such change need to have the Arrhenius activation energy (E) and the Arrhenius pre-exponential factor (A) defined.
 - During the heat simulation, you will see an illustration of the temperature in your modelled cube.
 - After this step, you will be shown three figures with the position of virtual temperature sensors in the cube, the temperature evolution at these positions, and whether there was some chemical change (based on the Arrhenius integral) in your cube.

5. (Optional) Calculate fluorescence light distribution (check out "Example4_FluorescenceAndImaging.m", section "%% Fluorescence Monte Carlo")
 - To be able to run this step, your geometry definitions needs to include the fluorescence wavelength, "wavelength_f".
 - "Example4_FluorescenceAndImaging.m" also contains the definitions for a "LightCollector" for both the incident light Monte Carlo and the fluorescence Monte Carlo. This enables simulating a light collection and imaging system, presented in a figure after this step. This is not required if you are only interested in the fluorescence light distribution in the cube.
 - After this step, you will be shown three figures with the fluorescence emitters distribution, the fluorescent light fluence rate and absorption within the cube. If you have chosen to use the LightCollector, you will additionally see an illustration of the geometry of your imaging system and the image of the fluorescent light.
  - You can also choose to show the light impinging on the LightCollector in a time-resolved manner. Check out "Example5_TimeTagging.m" on how to do so.

6. (Optional) Imaging: Model the intensity distribution received by a detector (check out ”Example4_FluorescenceAndImaging.m”, sections ”%% Monte Carlo simulation” and ”%% Fluorescence Monte Carlo”)
 - Include a LightCollector for the incident light Monte Carlo, the fluorescence Monte Carlo, or for both. The LightCollector can either be a pixel array detector or a single pixel detector (such as a fibre). You can also optionally make your detector time-resolved.

7. (Optional) Calculate the far field distribution of light escaping the simulation volume (see "Example6_FarField.m")

8. (Optional) Calculate the fluence rate matrix of only those photons which end up on the detector (see "Example7_CollectedFluenceRate.m")
 
9. (Optional) Programmatically assign values to the parameters
 - See "Example8_GeometryParametricSweep.m" and "Example9_MediaPropertyParametricSweep.m" on how to implement sweeps of the parameters.

## Contribution guidelines ##

If you want to report a bug or are missing a feature, report that as an issue on https://gitlab.gbar.dtu.dk/biophotonics/MCmatlab.

### Who do I talk to? ###

This repository is part of the DTU "biophotonics" team.
The main responsible is Anders Kragh Hansen: ankrh@fotonik.dtu.dk
