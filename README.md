# Readme and User Manual for MCmatlab

## What is MCmatlab?
MCmatlab is a Monte Carlo simulation for modeling light propagation in a 3D voxel space. Fluorescence can optionally be simulated after simulation of the excitation light. Included is also a finite element simulation for temperature increase and heat diffusion in the same voxel space.

Primarily targeted for tissue optics, but can be used in any environment with turbid media in which the wave nature of light (interference phenomena) is negligible and a ray-tracing model is appropriate.

You can find the latest version of MCmatlab on https://github.com/ankrh/MCmatlab. If you publish results obtained with this software, we would be thankful if you cited its accompanying article: doi:10.1117/1.JBO.23.12.121622

## License
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
- Windows 7, macOS 10.12 (Sierra) or newer or Linux
- MATLAB R2018a or newer
- (For GPU accelerated computation) A Windows PC with a CUDA-enabled graphics card and the MATLAB Parallel Computing Toolbox

### Helper files:
All the helper functions needed for running MCmatlab are located in the folder (MATLAB package) "+MCmatlab", and its parent folder therefore has to be on your MATLAB path when trying to run any MCmatlab model files. You will not need to modify any of the contained files. We recommend keeping the model files in a folder together with the "+MCmatlab" folder, such that you don't need to manually add "+MCmatlab" to your MATLAB path.

## How do I use MCmatlab?
### Compilation
The folders include all the executables necessary, so you don't need to compile anything. If, however, you want to change the routine in either the MCmatlab.c source code or the finiteElementHeatPropagator.c source code (both located in the folder "+MCmatlab/src"), you will need to recompile the respective mex-files. Check out those two source-files on how to do so. The source code is written in such a way that it can be compiled as either C or C++ code using either GCC, MSVC, clang or as CUDA source code using NVCC.

### Model files:
In MCmatlab, you set up your model in a single m-file. You can find some examples to get you started in the root folder of MCmatlab. A complete list of parameters that can be set is provided below.

### Media properties:
The optical properties of all media are defined in a function at the end of the model file. The examples contain some media you can modify or copy-paste to your own model file, or you define your own media from scratch in your own models. Make sure that each media in a single model file has its distinct "j"-number, as this is the number you will refer to when building your model in the geometry function (see below).

Each medium can be specified with up to 10 properties (fields in the struct):
- `name` - A char array containing the name of the medium, used for the visualization.
- `mua` - The absorption coefficient in units of 1/cm.
- `mus` - The scattering coefficient in units of 1/cm.
- `g` - The scattering anisotropy (the mean cosine to the scattering angle), to be used in the Henyey-Greenstein phase function. It has no units.
- `customPhaseFunc` - An alternative to specifying `g`. This is a char array that is a MATLAB-evaluateable expression that specifies the scattering phase function as a function of theta. This function does not need to be normalized.
- `PY` - The power yield of fluorescence, that is, watts of fluorescence emitted relative to watts of excitation light absorbed. This field is mutually exclusive with `QY`. It has no units.
- `QY` - The quantum yield of fluorescence, that is, photons of fluorescence emitted relative to excitation photons absorbed. This is related to the power yield through a factor of model.MC.lambda/model.FMC.lambda. This field is mutually exclusive with `PY`. It has no units.
- `VHC` - The volumetric heat capacity of the medium in units of J/(cm^3 K).
- `TC` - The thermal conductivity of the medium in units of W/(cm K).
- `E` - The Arrhenius activation energy in units of J/mol.
- `A` - The Arrhenius pre-exponential factor in units of 1/s.

For Monte Carlo simulations, `name`, `mua`, `mus` and `g` need be specified for all media.

For fluorescence Monte Carlo simulations, the `PY` or `QY` field must additionally be specified for those media that fluoresce.

For heat simulations, `VHC` and `TC` must be specified for all media.

For heat simulations with calculation of tissue damage, `E` and `A` must be specified for those media that can be damaged.

### Building the model
You build each model in a separate m-file. Each model requires the first two and optionally more of the following steps:
#### 1. Build geometry
- Modify or add to "mediaPropertiesFunc" at the end of the model file to include the definitions of the media you're interested in. You may optionally include thermal and/or fluorescence properties.
- In your model file, specify the side lengths Lx, Ly, Lz of your model cuboid and the resolutions Nx, Ny, Nz.
- The "GeomFunc" you define at the end of the model file simply returns a 3D array containing the media definition for each voxel in the model cuboid. The media are referred to by a number, corresponding to the "j" defined in "mediaPropertiesFunc" at the end of the model file.
- You may optionally import your geometry from STL files as shown in example 18.
- After this step, when calling "plot(model,'G')", you will be shown a figure with the geometry you defined.

#### 2. Calculate light distribution
- The section "%% Monte Carlo simulation" in the model file contains the definitions for the Monte Carlo simulation.
- After this step, when calling "plot(model,'MC')", you will be shown various figures; one with an overview of the optical properties, one with the absorbed light in the cuboid, one with the normalized fluence rate (NFR, roughly equivalent to normalized intensity or irradiance), and one showing the fluence rate at the cuboid boundaries (taking into account only the *exiting* photon packages, not the incident beam). The last figure will only be shown if some of your cuboid boundaries are "escaping" boundaries, i.e., boundary type 1 (all boundaries are escaping boundaries) or 2 (only the top is an escaping boundary).
- (Optional) If you set "model.MC.nExamplePaths" to an integer > 0, you will also be shown a plot with the paths of the number of photons you requested.

#### 3. (Optional) Include Fresnel reflection and refraction
- Check out "Example3_RefractionReflection.m" and "Example17_CurvedRefractionReflection.m" on how to implement changing refractive indices in your model.

#### 4. (Optional) Calculate fluorescence light distribution
- The section "%% Fluorescence Monte Carlo" contains the parameters for a fluorescence simulation, to be run after the distribution of excitation light has been calculated.
- To be able to run this step, your geometry definitions needs to include the fluorescence wavelength, "model.FMC.wavelength". Additionally, at least one of your media should be defined fluorescent, by defining its fluorescent yield Y in "mediaPropertiesFunc" at the end of your model file.
- After this step, when calling "plot(model,'FMC')", you will be shown the fluorescence emitter distribution, the fluorescent light fluence rate, the absorption within the cuboid, and the fluorescent light fluence rate at the cuboid boundaries.

#### 5. (Optional) Simulate heat distribution
- In the "%% Heat simulation" section, you put the parameters for simulation of the temperature distribution, to be run after the Monte Carlo step(s)
- To simulate the heat diffusion, all media involved in the simulation must have the thermal properties volumetric heat capacity (VHC) and thermal conductivity (TC) defined in "mediaPropertiesFunc" at the end of your model file.
- If you want to model chemical changes such as tissue damage, the media that might undergo such change need to additionally have the Arrhenius activation energy (E) and the Arrhenius pre-exponential factor (A) defined.
- During the heat simulation, you will see an illustration of the temperature in your modelled cuboid.
- After this step, when calling "plot(model,'HS')", you will be shown various figures with the position of virtual temperature sensors in the cuboid, the temperature evolution at these positions, and whether there was some chemical change (based on the Arrhenius integral) in your cuboid.
- You can run multiple heat simulations defining different successive pulse trains as shown in "Example11_MultipleHeatSims.m"

#### 6. (Optional) Imaging: Model the intensity distribution received by a detector/light collector
- Include a light collector for the incident light Monte Carlo, the fluorescence Monte Carlo, or for both. The light collector can either be a pixel array detector or a single pixel detector (such as a fibre). You can also optionally make your detector time-resolved.
- After running and calling either "plot(model,'MC')" or "plot(model,'FMC')", you will see an illustration of the geometry of your imaging system and the image of the scattered and/or fluorescence light.
- You can also choose to show the light impinging on the light collector in a time-resolved manner. Check out "Example6_TimeTagging.m" on how to do so.

#### 7. (Optional) Calculate the far field angular distribution of light escaping the simulation volume
- See "Example7_FarField.m". Note that the far field does not include "killed" photons, that is, photons that have exited at a boundary where the refractive index is different to 1, or photons that exit on boundaries that are not "escaping".

#### 8. (Optional) Calculate the normalized fluence rate array of only those photons which end up on the detector
- See "Example8_CollectedFluenceRate.m".

#### 9. (Optional) Programmatically assign values to the parameters
- See "Example9_GeometryParametricSweep.m" and "Example10_MediaPropertyParametricSweep.m" on how to implement sweeps of the parameters.

#### 10. (Optional) Define fluence rate-, temperature- or fractional damage-dependent optical or thermal properties
- See Examples 12-15.
- By specifying the optical or thermal properties as a string (technically a char array) defining the formula as a function of fluence rate (FR), temperature (T) or fractional damage (FD), MCmatlab will take the dependence into account with an iterative approach.

#### 11. (Optional) Enable CUDA GPU acceleration
- See Example16_CUDAacceleration.m.
- If you have an Nvidia graphics card with compute capability at least 3.0, you can set useGPU = true to enable running on the GPU, greatly speeding up both the Monte Carlo and heat simulations.

### List and explanation of input parameters
In the following we assume that the model object variable has been named "model". In principle, it could be given any name you want.
#### Geometry parameters
`model.G.silentMode`
[-]
(Default: False) If true, MCmatlab will not make any command window text outputs during the geometry initialization.

`model.G.nx`, `model.G.ny`, `model.G.nz`
[-]
The number of bins in the x, y, and z direction, respectively. Higher numbers mean a more finely resolved result, but will also require much longer runtimes to reach a low noise level.

`model.G.Lx`, `model.G.Ly`, `model.G.Lz`
[cm]
The side lengths of the simulation cuboid, measured in cm.

`model.G.mediaPropertiesFunc`
[-]
A MATLAB function handle. This should be set to refer to the function at the end of the model file in which the media properties are defined. In principle you could have several different media property functions in your model file and the one referred to in this property is the one that MCmatlab uses, but in practice most model files will be built using only one media property function.

`model.G.mediaPropParams`
[user-defined units]
A user-defined cell array that you can use to contain all sorts of inputs you would like to use inside the mediaPropertiesFunc. This cell array is useful when doing parametric sweeps (many MC simulations in a for or while loop) in which the media properties vary.

`model.G.geomFunc`
[-]
(Similar to mediaPropertiesFunc above)
A MATLAB function handle. This should be set to refer to the function at the end of the model file in which the model geometry is defined. In principle you could have several different geometry functions in your model file and the one referred to in this property is the one that MCmatlab uses, but in practice most model files will be built using only one geometry function.

`model.G.geomFuncParams`
[user-defined units]
(Similar to mediaPropParams above)
A user-defined cell array that you can use to contain all sorts of inputs you would like to use inside the geomFunc. This cell array is useful when doing parametric sweeps (many MC simulations in a for or while loop) in which the geometry varies.

#### (Excitation) Monte Carlo parameters
`model.MC.useGPU`
[-]
(Default: False)
If true, will run the MC simulation on GPU. If false, will use CPU.

`model.MC.GPUdevice`
[-]
(Default: 0)
The GPU device to use for GPU accelerated computations. 0 corresponds to the first device. For systems with multiple GPUs, you can launch simulations to different GPUs by varying this property and with the use of, e.g., a MATLAB parallel-for (parfor) loop.

`model.MC.simulationTimeRequested`
[minutes]
(Mutually exclusive with model.MC.nPhotonsRequested)
The time to run the MC simualtion for, in minutes. The number of photos launched will vary from run to run.

`model.MC.nPhotonsRequested`
[-]
(Mutually exclusive with model.MC.simulationTimeRequested)
The number of photons to launch. The execution time will vary from run to run.

`model.MC.silentMode`
[-]
(Default: False)
If true, MCmatlab will not make any command window text outputs during the MC simulation.

`model.MC.useAllCPUs`
[-]
(Default: False)
(Has no effect on MacOS or if model.MC.useGPU = true)
If true, MCmatlab will launch a number of threads equal to the number of CPU logical processors on the system. This will lead to the fastest execution but may make the system slow to respond to user inputs (mouse clicks etc.) for the duration of the simulation.
If false, MCmatlab will still be multithreaded but will leave one CPU logical processor unused. This will yield almost the same speed of execution but will make the system more responsive while the simulation is running.

`model.MC.calcNFR`
[-]
(Default: True)
If true, will calculate and store the normalized fluence rate (NFR) 3D array. The array takes 8\*nx\*ny\*nz bytes of memory (and disk space, if the resulting model file is saved to disk subsequently)
For some simulations in which you are only interested, for example, in the detector image/power, you might not need the NFR and could therefore set this to false.

`model.MC.calcNFRdet`
[-]
(Default: False)
If true, will additionally calculate a normalized fluence rate 3D array similar to the ordinary NFR array, but only counting the areas passed through by those photons that end up registered by the detector/light collector.

`model.MC.nExamplePaths`
[-]
(Default: 0)
If set to a positive integer N, will store the paths of the first N photons for subsequent visualization during the plotting. This is useful to visualize how the photons move around the cuboid.

`model.MC.farFieldRes`
[pixels]
(Default: 0)
If set to a positive integer N, MCmatlab will calculate the far field direction angles of the escaping photons and store them in a 2D NxN array (theta,phi on the unit sphere) for later visualization.

`model.MC.matchedInterfaces`
[-]
(Default: True)
If true, will set all refractive indices to 1, disregarding any refractive indices specified in the media properties function. In this case, no Fresnel reflection and refraction will be simulated.

`model.MC.smoothingLengthScale`
[cm]
(Only used if matchedInterfaces = false)
The characteristic length scale to use for smoothing of the Sobel filter applied to the normal vectors of the non-matched interfaces within the geometry. See
https://www.spiedigitallibrary.org/journals/journal-of-biomedical-optics/volume-25/issue-02/025001/Modeling-voxel-based-Monte-Carlo-light-transport-with-curved-and/10.1117/1.JBO.25.2.025001.full
for a description of the technique used.
The use of such smoothing enables simulation of reflection and refraction even on oblique interfaces, despite the geometry being defined on a rectangular cuboid mesh. See Example17_CurvedRefractionReflection.m.

`model.MC.boundaryType`
[-]
0: No escaping boundaries
1: All cuboid boundaries are escaping
2: Top cuboid boundary only is escaping
3: Top and bottom boundaries are escaping, while the side boundaries are cyclic
Photons that hit a boundary that is "escaping" and where the refractive index is equal to 1 will be considered an "escaped" photon, and may, depending on the position and direction, be registered by the detector/light collector.
When a photon encounters a boundary that is not escaping, it is allowed to propagate beyond the edge up to a distance of 2?? times the value of (Lx or Ly or Lz) before the photon is finally killed and no longer simulated.
Allowing photons to propagate beyond the cuboid boundaries enable the photons to potentially scatter back into the simulation volume.
The more boundaries are escaping, the faster the simulation will run because fewer photons outside the cuboid have to be simulated.
When a photon hits a cyclic boundary, it will exit and immediately enter the cuboid again on the opposite boundary, as if the cuboid is periodic. In other words, if your geometry and light source is periodically repeating in x and y (such as a simple layered model illuminated by an infinite plane wave), you need only simulate just one unit cell of this pattern and use cyclic x and y boundaries.

`model.MC.wavelength`
[nm]
The wavelength of the light. This is used only to be passed on as an input to the mediaPropertiesFunc because most (if not all) optical parameters are dependent on the wavelength of the light. The wavelengths does not actually change anything in the MC simulation itself aside from the change in mua, mus and g of the involved media. The only exception is when simulation fluorescence (see below) and specifying a quantum yield fluorescer.

`model.MC.lightSource.sourceType`
[-]
(Note that lightSource can be abbreviated LS in your code)
0: Pencil beam
A pencil beam is a beam with zero width and zero divergence.
1: Isotropically emitting line or point source
This option is the only one in which the light will originate from within the cuboid. All other options have photons packets starting at the top surface, z = 0.
2: Infinite plane wave
A beam with infinite width, but zero divergence.
3: Laguerre-Gaussian LG01 beam
A ray-based emulation of the donut-shaped Laguerre-Gaussian LG01 beam.
4: Radial-factorizable beam
A beam in which the focus beam profile can be written as I(x,y) = Ir\(r\), where r = sqrt(x\^2+y\^2), and similarly for the angular distribution. Circular Gaussian beams can be written in this form.
5: X/Y factorizable beam
A beam in which the focus beam profile can be written as I(x,y) = Ix(x)*Iy(y), and similarly for the angular distribution. Rectangular LED emitters can be written in this way, as well as Gaussian (circular and elliptical) beams.

`model.MC.lightSource.emitterLength`
[cm]
(Only used for model.MC.lightSource.sourceType = 1. Default: 0)
The length of the line emitter. If zero, the emitter is a point source.

`model.MC.lightSource.focalPlaneIntensityDistribution.radialDistr`
[-]
(Note that focalPlaneIntensityDistribution can be abbreviated FPID in your code)
(Only used for sourceType = 4)
Radial intensity distribution in the focal plane.
0: Top-hat
1: Gaussian
Array: Custom. Doesn't need to be normalized.

`model.MC.lightSource.focalPlaneIntensityDistribution.radialWidth`
[cm]
(Note that focalPlaneIntensityDistribution can be abbreviated FPID in your code)
(Only used for sourceType = 4)
Radial width of the intensity distribution in the focal plane.
If FPID.radialDistr is set to 0 or 1, this is the 1/e^2 radius
If FPID.radialDistr is set to an array, this is the half-width of the full distribution

`model.MC.lightSource.focalPlaneIntensityDistribution.XDistr`, `model.MC.lightSource.focalPlaneIntensityDistribution.YDistr`
[-]
(Note that focalPlaneIntensityDistribution can be abbreviated FPID in your code)
(Only used for sourceType = 5)
X and Y distribution of the intensity in the focal plane.
0: Top-hat
1: Gaussian
Array: Custom. Doesn't need to be normalized.

`model.MC.lightSource.focalPlaneIntensityDistribution.XWidth`, `model.MC.lightSource.focalPlaneIntensityDistribution.YWidth`
[cm]
(Note that focalPlaneIntensityDistribution can be abbreviated FPID in your code)
(Only used for sourceType = 5)
X and Y width of the intensity distribution in the focal plane.
If FPID.XDistr/FPID.YDistr is set to 0 or 1, this is the 1/e^2 radius
If FPID.XDistr/FPID.YDistr is set to an array, this is the half-width of the full distribution

`model.MC.lightSource.angularIntensityDistribution.radialDistr`
[-]
(Note that angularIntensityDistribution can be abbreviated AID in your code)
(Only used for sourceType = 4)
Radial distribution of the angular intensity.
0: Top-hat
1: Gaussian
2: Cosine (Lambertian)
Array: Custom. Doesn't need to be normalized.

`model.MC.lightSource.angularIntensityDistribution.radialWidth`
[rad]
(Note that angularIntensityDistribution can be abbreviated AID in your code)
(Only used for sourceType = 4 and if AID.radialDistr is not 2)
Radial width of the angular intensity distribution.
If AID.radialDistr is set to 0 or 1, this is the 1/e^2 half-angle.
If AID.radialDistr is set to an array, this is the half-angle of the full distribution.
For a diffraction limited Gaussian beam, this should be set to model.MC.wavelength\*1e-9/(pi\*model.MC.lightSource.focalPlaneIntensityDistribution.radialWidth\*1e-2))

`model.MC.lightSource.angularIntensityDistribution.XDistr`, `model.MC.lightSource.angularIntensityDistribution.YDistr`
[-]
(Note that angularIntensityDistribution can be abbreviated AID in your code)
(Only used for sourceType = 5)
X and Y distribution of the angular intensity.
0: Top-hat
1: Gaussian
2: Cosine (Lambertian)
Array: Custom. Doesn't need to be normalized.

`model.MC.lightSource.angularIntensityDistribution.XWidth`, `model.MC.lightSource.angularIntensityDistribution.YWidth`
[rad]
(Note that angularIntensityDistribution can be abbreviated AID in your code)
(Only used for sourceType = 5 and if AID.XDistr or AID.YDistr is not 2)
Radial width of the angular intensity distribution.
If AID.XDistr/AID.YDistr is set to 0 or 1, this is the 1/e^2 half-angle.
If AID.XDistr/AID.YDistr is set to an array, this is the half-angle of the full distribution.

`model.MC.lightSource.xFocus`, `model.MC.lightSource.yFocus`, `model.MC.lightSource.zFocus`
[cm]
The focus position x,y,z coordinates.

`model.MC.lightSource.psi`
[rad]
(Only used for sourceType = 5. Default: 0)
Axial rotation angle of the beam.

`model.MC.lightSource.theta`, `model.MC.lightSource.phi`
[rad]
Polar and azimuthal angle of the beam center axis.

`model.MC.depositionCriteria.minScatterings`, `model.MC.depositionCriteria.maxScatterings`, `model.MC.depositionCriteria.minRefractions`, `model.MC.depositionCriteria.maxRefractions`, `model.MC.depositionCriteria.minReflections`, `model.MC.depositionCriteria.maxReflections`, `model.MC.depositionCriteria.minInterfaceTransitions`, `model.MC.depositionCriteria.maxInterfaceTransitions`
[-]
The minimum and maximum number of scattering, refraction, reflection and interface transition events the photon packet must have undergone in order to deposit weight into the output arrays (normalized fluence rate, normalized boundary irradiances, image and far field). There is an implied AND between all the criteria. Note that the underlying absorption and propagation of all photon packets are independent of the specified deposition criteria.

`model.MC.P`
[W]
(Only used for the thermal simulations and for any Monte Carlo simulations in which the optical or thermal parameters depend on the fluence rate (FR), fractional damage (FD) or temperature (T))
Incident pulse peak power. In case of infinite plane waves, this is only the power incident upon the cuboid's top surface.

`model.MC.FRinitial`
[W/cm^2]
(Only used for MC simulations in which the optical properties depend on the fluence rate (FR))
Initial guess for the fluence rate distribution in, to be used for iterative fluence rate dependent simulations.

`model.MC.FRdepIterations`
[-]
(Only used for MC simulations in which the optical properties depend on the fluence rate (FR))
The number of iterations (times the MC simulation needs to be run) each time a new calculation of the MC outputs is requested.
Test for yourself how high this number needs to be for your model. In my tests, 20 was enough.

`model.MC.useLightCollector`
[-]
(Default: False)
A logical (boolean) to specify whether the MC simulation should keep track of whether photons hit the light collector/detector.


(The following LC properties are only used if useLightCollector = true)

`model.MC.lightCollector.f`
[cm]
(Note that lightCollector can be abbreviated LC in your code)
If the light collector is supposed to be an objective lens, this is its focal length.
If the light collector is supposed to be a fiber tip this should be set to Inf (infinity)

`model.MC.lightCollector.x`, `model.MC.lightCollector.y`, `model.MC.lightCollector.z`
[cm]
(Note that lightCollector can be abbreviated LC in your code)
If model.MC.lightCollector.f is finite (light collector is an objective lens), this is the x, y, z position of the center of the focal plane of the objective lens.
If model.MC.lightCollector.f is infinite (light collector is a fiber tip), this is the x, y, z position of the center of the fiber tip.

`model.MC.lightCollector.theta`, `model.MC.lightCollector.phi`
[rad]
(Note that lightCollector can be abbreviated LC in your code)
Polar and azimuthal angle of direction that the light collector is facing

`model.MC.lightCollector.NA`
[-]
(Note that lightCollector can be abbreviated LC in your code)
(Only used for infinite f)
Fiber numerical aperture.

`model.MC.lightCollector.diam`
[cm]
(Note that lightCollector can be abbreviated LC in your code)
Diameter of the light collector aperture. For an ideal thin lens with numerical aperture NA, this is 2f\*tan(asin(NA)).

`model.MC.lightCollector.fieldSize`
[cm]
(Note that lightCollector can be abbreviated LC in your code)
(Only used for finite f)
Field size of the imaging system, that is, the diameter of area in object plane that gets imaged.

`model.MC.lightCollector.res`
[pixels]
(Note that lightCollector can be abbreviated LC in your code)
(Only used for finite f)
X and Y resolution of light collector, only used for finite f

`model.MC.lightCollector.tStart`, `model.MC.lightCollector.tEnd`
[s]
(Note that lightCollector can be abbreviated LC in your code)
Start and end time of the detection time-of-flight interval. Note that photons arriving before tStart will still be registered but all binned together in a single bin containing all "too-early" photons. Likewise, all photons arriving after tEnd are binned into a "too-late" bin.

`model.MC.lightCollector.nTimeBins`
[-]
(Note that lightCollector can be abbreviated LC in your code)
(Default: 0)
Number of bins between tStart and tEnd. If zero, the measurement is not time-resolved.

#### Fluorescence Monte Carlo parameters
The following properties exist for fluorescence Monte Carlo simulations, and they work the same as for regular MC simulations:
`FMC.useGPU`, `FMC.GPUdevice`, `FMC.simulationTimeRequested`, `FMC.nPhotonsRequested`, `FMC.silentMode`, `FMC.useAllCPUs`, `FMC.calcNormalizedFluenceRate`, `FMC.calcNormalizedFluenceRate_detected`, `FMC.nExamplePaths`, `FMC.farFieldRes`, `FMC.matchedInterfaces`, `FMC.smoothingLengthScale`, `FMC.boundaryType`, `FMC.wavelength`, `FMC.useLightCollector`, `FMC.lightCollector.x`, `FMC.lightCollector.y`, `FMC.lightCollector.z`, `FMC.lightCollector.theta`, `FMC.lightCollector.phi`, `FMC.lightCollector.f`, `FMC.lightCollector.diam`, `FMC.lightCollector.fieldSize`, `FMC.lightCollector.NA`, `FMC.lightCollector.res`

#### Heat solver parameters
`model.HS.useGPU`
[-]
(Default: false)
Use CUDA acceleration for the heat simulation solver for NVIDIA GPUs

`model.HS.silentMode`
[-]
(Default: false)
If true, MCmatlab will not make any command window text outputs during the heat simulation.

`model.HS.useAllCPUs`
[-]
(Default: false)
(Has no effect on MacOS or if model.MC.useGPU = true)
If true, MCmatlab will launch a number of threads equal to the number of CPU logical processors on the system. This will lead to the fastest execution but may make the system slow to respond to user inputs (mouse clicks etc.) for the duration of the simulation.
If false, MCmatlab will launch will still be multithreaded but will leave one CPU logical processor unused. This will yield almost the same speed of execution but will make the system more responsive while the simulation is running.

`model.HS.makeMovie`
[-]
(Default: false)
(Only used if silentMode = false)
If true, MCmatlab will record each frame of the volumetric temperature plot and get ready to save the frames to a movie file.

`model.HS.deferMovieWrite`
[-]
(Default: false)
(Only used if silentMode = false and makeMovie = true)
If false, MCmatlab will save the recorded frames to a movie file. If true, the frames are still stored within the model object so that subsequent heat simulations can be concatenated to them before saving.

`model.HS.largeTimeSteps`
[-]
(Default: false)
If true, calculations will be faster, but some voxel temperatures may be slightly less precise. Test for yourself whether this precision is acceptable for your application.

`model.HS.heatBoundaryType`
[-]
The type of boundary condition to use for the heat solver.
0: Insulating boundaries
1: Constant-temperature boundaries (heat-sinked)

`model.HS.nPulses`
[-]
Number of consecutive pulses, each with an illumination phase and a diffusion phase. If simulating only illumination or only diffusion, use nPulses = 1.

`model.HS.durationOn`, `model.HS.durationOff`
[s]
Pulse on-duration and off-duration

`model.HS.durationEnd`
[s]
Non-illuminated relaxation time to add to the end of the simulation to let temperature diffuse after the pulse train.

`model.HS.Tinitial`
[deg C]
Initial temperature, can be scalar or 3D array of dimensions nx, ny, nz.

`model.HS.plotTempLimits`
[deg C]
An array of two values designating the expected range of temperatures, used only for setting the color scale in the plot.

`model.HS.nUpdates`
[-]
Number of times data is extracted for plots during each pulse. A minimum of 1 update is performed in each phase (2 for each pulse consisting of an illumination phase and a diffusion phase)

`model.HS.mediaPropRecalcPeriod`
[-]
(Used only if optical or thermal parameters have been specified to be dependent on temperature (T))
Every time this many updates have passed, the media properties will be recalculated (including, if needed, re-running MC and FMC steps). Check that this is low enough to ensure convergence of your results.

`model.HS.slicePositions`
[-]
(Default: [0.5 1 1])
An array of three values, describing the fractional x, y, z positions on a scale from 0 to 1 to place the volumetric slices at.

`model.HS.tempSensorPositions`
[cm]
(Default: empty array)
A 2D array with three columns.
Each row describes the x, y, z coordinates of a temperature sensor placed within the cuboid, at which the temperature will be recorded for subsequent plotting as a function of time.
Leave the matrix empty ([]) to disable temperature sensors.

### List and explanation of output properties
As in the previous section, we assume that the model object variable has been named "model". In principle, it could be given any name you want.

#### Geometry properties
`model.G.dx`, `model.G.dy`, `model.G.dz`
[cm]
The side lengths of the voxels. They are simply calculated as the side length of the simulation cuboid divided by the number of bins, for example, model.G.Lx/model.G.nx.

`model.G.x`, `model.G.y`, `model.G.z`
[cm]
1D arrays with the x, y or z coordinate values corresponding to the centers of the voxels
You can use these for easy plotting of slices, see for example the entry on NFR in the MC section below.

`model.G.M_raw`
[-]
3D (xyz) array in coordinate order of the tissue indices as calculated from the geomFunc. This is the array that is plotted in the figure called Geometry Illustration.

#### (Excitation) Monte Carlo properties
`model.MC.simulationTime`
[min]
The actual simulation time of the most recent Monte Carlo simulation run. If you did not define model.MC.nPhotonsRequested, then the simulation is timed and model.MC.simulationTime should be very close to `model.MC.simulationTimeRequested.

`model.MC.nPhotons`
[-]
The actual number of photon packets launched in the most recent Monte Carlo simulation run. If you defined model.MC.nPhotonsRequested, then model.MC.nPhotons will be equal to that number.

`model.MC.mediaProperties`
[-]
A struct that contains the media properties as evaluated at the specified excitation Monte Carlo wavelength, model.MC.wavelength.

`model.MC.M`
[-]
A 3D (xyz) array. If the media properties do not contain dependences on fluence rate (FR), temperature (T) or fractional damage (FD), then this array will simply be equal to model.G.M_raw. Otherwise, model.MC.M will be the "split" version of model.G.M_raw, in which those media that have dependences have been split into the specified number of sub-media with different optical/thermal properties for the purposes of the simulations.

`model.MC.interfaceNormals`
[-]
A 2D array of dimension (2 x (nx*ny*nz)) where local interface normal vector angles are stored. Each column contains the theta, phi angles for one voxel's interface normal vector. Used for calculated of refraction and reflection on oblique and curved interfaces between media with different refractive index.

`model.MC.examplePaths`
[cm] and [-]
A 2D array which stores the example paths' coordinates and remaining weight. Each column is a data point where the first three rows are the x, y, z position of the photon packet and the fourth row is the remaining weight (energy) of the photon packet. Different photon trajectories are separated by a column of NaNs.

`model.MC.normalizedFluenceRate`
[W/cm^2/W.incident]
(Note that normalizedFluenceRate can be abbreviated NFR in your code)
A 3D (xyz) array of normalized fluence rate values. This is what's plotted in figure 4. If you want to plot a slice out of this distribution, you could for example do it with these three lines of code:
- `NFRslice = squeeze(model.MC.NFR(51,:,:));` Extracts the slice corresponding to the 51st x value. `squeeze` is to remove the x singleton dimension, giving you a 2D (yz) array of normalized fluence rates.
- `figure(100);`
- `imagesc(model.G.y,model.G.z,NFRslice.');` We can use the 1D coordinate arrays from model.G for the y and z coordinate values. The `.'` is there because Mathworks has designed the imagesc function to plot the first coordinate along the vertical axis and the second coordinate along the horizontal axis, which is the opposite of what we want.

`model.MC.normalizedFluenceRate_detected`
[W/cm^2/W.incident]
(Note that normalizedFluenceRate can be abbreviated NFR in your code)
A 3D (xyz) array of normalized fluence rate values like model.MC.NFR, but only counting those photons that ended up on the light collector.

`model.MC.farField`
[W/sr/W.incident]
A 2D (theta,phi) array of normalized radiant intensity values for the light that escaped the simulation cuboid. The axes are the polar and azimuthal angles, described below.

`model.MC.farFieldTheta`
[rad]
A 1D array of polar angles corresponding to the first dimension of model.MC.farField.

`model.MC.farFieldPhi`
[rad]
A 1D array of azimuthal angles corresponding to the second dimension of model.MC.farField.

`model.MC.normalizedIrradiance_xpos`, `model.MC.normalizedIrradiance_xneg`, `model.MC.normalizedIrradiance_ypos`, `model.MC.normalizedIrradiance_yneg`, `model.MC.normalizedIrradiance_zpos`, `model.MC.normalizedIrradiance_zneg`
[W/cm^2/W.incident]
(Note that normalizedIrradiance can be abbreviated NI in your code)
2D (xy, xz or yz) arrays of normalized irradiances of the light hitting the boundaries (either escaping or killed). For example, the property named `model.MC.NI_xpos` designates the light that hits the surface that is orthogonal to the x axis and placed on the positive x side, whereas the `model.MC.NI_xneg` refers to the surface orthogonal to the x axis and placed on the negative x side. model.MC.boundaryType determine which of these arrays are calculated.
`model.MC.NI_zneg` is of special interest to users interested in calculating the reflectance, which can be found as the integral over the array:
- `R = model.G.dx*model.G.dy*sum(model.MC.NI_zneg(:));`

`model.MC.lightCollector.image`
[W/cm^2/W.incident]
If `model.MC.lightCollector.res == 1`, this is a scalar with the normalized power registered on the light collector.
If `model.MC.lightCollector.res > 1`, this is a 2D (XY) array of normalized irradiances registered on the light collector. This array describes the light distribution you would get on a camera looking at the cuboid through an objective lens.

#### Fluorescence Monte Carlo properties
All the output properties described above for excitation Monte Carlo are also provided for fluorescence Monte Carlo in the model.FMC object.

#### Heat solver properties
`model.HS.M`
[-]
A 3D (xyz) array. If the media properties do not contain dependences on fluence rate (FR), temperature (T) or fractional damage (FD), then this array will simply be equal to model.G.M_raw. Otherwise, model.HS.M will be the "split" version of model.G.M_raw, in which those media that have dependences have been split into the specified number of sub-media with different optical/thermal properties for the purposes of the simulations.

`model.HS.T`
[deg C]
A 3D (xyz) array with the final temperature values of the simulation.

`model.HS.Omega`
[-]
A 3D (xyz) array with the final Omega values of the simulation. Omega describes the thermally induced tissue damage calculated according to the Arrhenius damage integral.

`model.HS.maxMediaTemps`
[deg C]
A 1D array that contains the maximal temperatures obtained in each tissue type involved in the simulation.

`model.HS.sensorsTimeVector`
[s]
A 1D array of time values that should be used for plotting the sensorTemps data (see below).

`model.HS.sensorTemps`
[deg C]
A 2D array of temperatures measured at the specified sensor positions. Each row corresponds to one sensor. Use the model.HS.sensorsTimeVector as the time values when plotting this data.

## Contribution guidelines
You are very welcome to clone the repository and add features of your own. If you want to report a bug or are missing a feature, shoot me an email, see below.

## Who do I talk to?
The main person responsible for MCmatlab is Anders Kragh Hansen: ankrh@fotonik.dtu.dk
