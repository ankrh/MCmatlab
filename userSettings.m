function userSettings

% In this code, user settings can be entered for the input spectrum to be
% simulated. The user defined variables is: dnm, nm_min, nm_max, dF, F_min,
% F_max and filename. Further, the text file where the spectrum is stored,
% defined by filename, is loaded and saved to the parameters nm and Power.

%Running the script will save these parameters to the file "Input_spectrum.mat", 
%and overwrite the previous settings saved  in that file.

%If any of the user settings are typed in incorrectly, fx. being out of range, 
%a message will occur in the command window stating the error. In this case, 
%the settings will not be saved. The wrong setting should be corrected and 
%the script run again.


%% User setting: Define the interval of wavelengths the input spectrum covers.

dnm = 0*5; %Wavelength step size of the spectrum. Needs to be a multiple of 5. [nm];
nm_min = 850; %Minimum wavelength of the filter. No lower than 300+dnm [nm];
nm_max = 850; %Maximum wavelength of the filter. No higher than 1000-dnm [nm];


%% User setting: Define the pulse. 

%Set the fluence of the pulse. 
%The selected max. and min. fluences, and the step size should be chosen so
%that it covers no more than 10 data points, for the sake of computation
%time. If only one value should be simulated, F_min and F_max should be put
%to the same number.

dF = 0.2; %Fluence step size. [J/cm2]
F_min = 20; %Minimum fluence to be simulated.
F_max = 20; %Maximum fluence to be simulated.

%Set the pulse duration in seconds
pulse = 10e-3; %pulse duration. Default 10 ms. [s]


%% User setting: Define the input spectrum file

%name of the spectrum file. Extension (.txt) should not be written as part
%of the filename. File must be present in current folder or a subfolder.
filename = 'VL+_norm0'; 
directoryPath = 'C:\Users\Kira Schmidt\Desktop\mcxyz';

%% Do not edit this section. (Calling the spectrum)
%Checks if the file specified in filename is in the folder or subfolder

AllFiles = dirrec(char(fileparts(pwd)),[filename '.txt*']); %file directory

%Checks that the file exists, and is unique
if length(AllFiles)>1
    disp('ERROR: The name specified for the file to load the spectrum from was not unique')
    return
elseif isempty(AllFiles) == 1
    disp('ERROR: The name specified for the file to load the spectrum from was not found')
    return
end

%Reads the file
fid = fopen(char(AllFiles),'r');
A = textscan(fid, '%f %f', 'HeaderLines', 0,'Delimiter','\t');
fclose(fid);

nm = A{1};
Power = A{2};

clear AllFiles fid A

%% Do not edit this section. (Error checking and file saving)

%Wavelength error checking
if (dnm/5) ~= round(dnm/5)
    disp('ERROR: Wavelength step size, dnm, is not an integer of 5')
    return
elseif nm_min > nm_max
    disp('ERROR: Minimum wavelength is bigger than maximum wavelength')
    return
elseif dnm < 5 || dnm > (nm_max - nm_min)
    disp('ERROR: Wavelength step size, dnm, out of range')
    return
elseif nm_min < 300 + dnm
    disp('ERROR: Minimum wavelength is too low')
    return
elseif nm_max > 1000 - dnm
    disp('ERROR: Maximum wavelength is too high')
    return
end

%Fluence error checking
if dF <= 0
    disp('ERROR: Fluence step size, dF, out of range')
    return
elseif F_min > F_max
    disp('ERROR: Minimum fluence bigger than maximum fluence')
    return
elseif dF >= (F_max - F_min) && (F_max - F_min) > 0
    disp('ERROR: Fluence step size, dF, out of range')
    return
elseif floor((F_max - F_min)/dF + 1) > 10
    disp('ERROR: Too many fluence data points')
    return
end

%Checks if all appropriate values are numbers
if isnumeric(dnm) == 0 || isnumeric(nm_min) == 0 || isnumeric(nm_max) == 0
    disp('ERROR: An entered wavelength setting is not a number')
    return
elseif isnumeric(dF) == 0 || isnumeric(F_min) == 0 || isnumeric(F_max) == 0
    disp('ERROR: An entered fleunce setting is not a number')
    return
end

dnm = 2*dnm;

save([directoryPath 'Input_spectrum'])
