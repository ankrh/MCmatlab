# MCmatlab

## What is this repository for?

A Monte Carlo simulation for modeling light propagation in a 3D voxel space.
Fluorescence can optionally be simulated after simulation of the excitation light.
Included is also a finite element simulation for temperature increase and heat diffusion in the same voxel space.
Primarily targeted for tissue optics, but can be used in any environment.

You can find the latest version of MCmatlab and more information in the wiki on https://gitlab.gbar.dtu.dk/biophotonics/MCmatlab.
If you publish results obtained with this software, we would be thankful if you cited its accompanying article: doi:10.1117/1.JBO.23.12.121622

## LICENSE

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

## How do I get set up?

### Requirements:
- Windows 7 or macOS 10.12 (Sierra) or newer
- MATLAB R2017a or newer

### Helper files:
All the helper functions needed for running MCmatlab are located in the folder "helperfuncs", which therefore has to be on your MATLAB path. You will not need to modify any of these files. The example model files automatically add "helperfuncs" to the MATLAB path as the first command, and you should keep that practice also in your own model files.

### Model files:
In MCmatlab, you set up your model in a single m-file. You can find a few examples to get you started in the root folder of MCmatlab. Once you're familiar with those, you can check the model file "Template.m" for all the possible switches you might use, but it is itself not a valid model file, as it also contains some mutually exclusive switches.

### Media properties:
The optical properties of all media are defined in the end of the model file. The examples contain some media you can modify or copy-paste to your own model file, or you define your own media from scratch in your own models. Make sure that each media in a single model file has its distinct "j"-number, as this is the number you will refer to when building your model.

### Versions of MCmatlab:
There currently (February 2020) exist two versions of MCmatlab, the old version (R-version), where the media-properties are specified in a file separate from the model file, and the new version (S-version), where the media properties are stored in the model file itself, which allows more flexible beam definitions, and which can use CUDA if you have a CUDA-enabled graphics card installed. Check in any of the example files whether you have a function "function mediaProperties = mediaPropertiesFunc(wavelength,parameters)". If you have this, you are on the new S-version. If not, you are on the old R-version and might want to download the new S-version from https://gitlab.gbar.dtu.dk/biophotonics/MCmatlab. You will then need to upgrade your old model files, which is straightforward. Read "Changes to model file commands, January 2020.txt" on how to do so. (If you prefer the old R-version of MCmatlab to run your old model files, you can download it from https://gitlab.gbar.dtu.dk/biophotonics/MCmatlab/repository/R1908b/archive.zip)

## How do I use MCmatlab?
### Compilation

The folders include all the executables necessary, so you don't need to compile anything. If, however, you want to change the routine in either the MCmatlab.c source code or the finiteElementHeatPropagator.c source code (both located in the folder "src"), you will need to recompile the respective mex-files. Check out those two source-files on how to do so. The source code is written in such a way that it can be compiled as either C or C++ code using either GCC or MSVC, or as CUDA source code using NVCC.

### Building the model
You build a model in a separate m-file. Each model requires the first two and optionally more of the following steps:

1. Build geometry (check out "Example1_StandardTissue.m", section "%% Geometry definition")
 - Modify or add to "mediaPropertiesFunc" at the end of the model file to include the definitions of the media you're interested in. You may optionally include thermal and/or fluorescence properties.
 - In your model file, specify the geometry of your model.
 - The "GeomFunc" you define in the model file simply defines a 3D-matrix containing the media definition for each voxel in the model cube. The media are referred to by a number, corresponding to the "j" defined in "mediaPropertiesFunc" at the end of the model file.
 - After this step, you will be shown a figure with the geometry you defined.

2. Calculate light distribution (check out "Example1_StandardTissue.m", section "%% Monte Carlo simulation")
 - This section in the model file contains the definitions for the Monte Carlo simulation.
 - After this step, you will be shown four figures; one with an overview of the optical properties, one with the normalized fluence rate, one with the absorbed light in the cube, and one showing the fluence rate at the cuboid boundaries (taking into account only the *exiting* photon packages, not the incident beam). The last figure will only be shown if some of your cuboid boundaries are "escaping" boundaries, i.e., boundary type 1 (all boundaries are escaping boundaries) or 2 (top boundary only is escaping boundary).

3. (Optional) Include Fresnel reflection and refraction
 - Check out "Example3_RefractionReflection.m" on how to implement changing refractive indices in your model.

4. (Optional) Simulate heat distribution (check out "Example4_BloodVessel.m", section "%% Heat simulation")
 - To simulate the heat diffusion, all media involved in the simulation must have the thermal properties volumetric heat capacity (VHC) and thermal conductivity (TC) defined in "mediaPropertiesFunc" at the end of your model file.
 - If you want to model chemical changes such as tissue damage, the media that might undergo such change need to have the Arrhenius activation energy (E) and the Arrhenius pre-exponential factor (A) defined.
 - During the heat simulation, you will see an illustration of the temperature in your modelled cube.
 - After this step, you will be shown three figures with the position of virtual temperature sensors in the cube, the temperature evolution at these positions, and whether there was some chemical change (based on the Arrhenius integral) in your cube.
 - You can run multiple heat simulations defining different successive pulse trains as shown in "Example11_MultipleHeatSims.m"

5. (Optional) Calculate fluorescence light distribution (check out "Example5_FluorescenceAndImaging.m", section "%% Fluorescence Monte Carlo")
 - To be able to run this step, your geometry definitions needs to include the fluorescence wavelength, "model.FMC.wavelength". Additionally, at least one of your media should be defined fluoerscent, by defining its fluorescent yield Y in "mediaPropertiesFunc" at the end of your model file.
 - After this step, you will be shown the fluorescence emitter distribution, the fluorescent light fluence rate, the absorption within the cube, and the fluorescent light fluence rate at the cuboid boundaries.

6. (Optional) Imaging: Model the intensity distribution received by a detector/light collector (check out ”Example5_FluorescenceAndImaging.m”, sections ”%% Monte Carlo simulation” and ”%% Fluorescence Monte Carlo”)
 - Include a light collector for the incident light Monte Carlo, the fluorescence Monte Carlo, or for both. The light collector can either be a pixel array detector or a single pixel detector (such as a fibre). You can also optionally make your detector time-resolved.
 - After running, you will see an illustration of the geometry of your imaging system and the image of the incident and/or fluorescence light.
 - You can also choose to show the light impinging on the light collector in a time-resolved manner. Check out "Example6_TimeTagging.m" on how to do so.

7. (Optional) Calculate the far field distribution of light escaping the simulation volume (see "Example7_FarField.m")

8. (Optional) Calculate the fluence rate matrix of only those photons which end up on the detector (see "Example8_CollectedFluenceRate.m")

9. (Optional) Programmatically assign values to the parameters
 - See "Example9_GeometryParametricSweep.m" and "Example10_MediaPropertyParametricSweep.m" on how to implement sweeps of the parameters.

10. (Optional) Define fluence rate-, temperature- or fractional damage-dependent optical or thermal properties (See Examples 12-15)
 - By specifying the optical or thermal properties as a string (technically a char array) defining the formula as a function of fluence rate (FR), temperature (T) or fractional damage (FD), MCmatlab will take the dependence into account with an iterative approach.

11. (Optional) Enable CUDA GPU acceleration (See Example16_CUDAacceleration.m)
 - If you have an Nvidia graphics card with compute capability at least 3.0, you can set useGPU = true to enable running on the GPU, greatly speeding up both the Monte Carlo and heat simulations.


## Contribution guidelines

If you want to report a bug or are missing a feature, report that as an issue on https://gitlab.gbar.dtu.dk/biophotonics/MCmatlab.

### Who do I talk to?

This repository is part of the DTU "biophotonics" team.
The main responsible is Anders Kragh Hansen: ankrh@fotonik.dtu.dk
