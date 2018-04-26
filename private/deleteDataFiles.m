function overwrite = deleteDataFiles(name)

fileNames = {...
    ['./Data/' name '.mat'] ...
    ['./Data/' name '_MCoutput.mat'] ...
    ['./Data/' name '_MCoutput_fluorescence.mat'] ...
    ['./Data/' name '_heatSimoutput.mat'] ...
    ['./Data/' name '_heatSimoutput.mp4'] ...
    };

caller = dbstack(1);

switch caller.name
    case 'makeTissue'
        fileNames = fileNames;
    case 'runMonteCarlo'
        fileNames = fileNames(2:end);
    case 'runMonteCarloFluorescence'
        fileNames = fileNames(3);
    case 'simulateHeatDistribution'
        fileNames = fileNames(4:end);
end
    

for fileNameIndex = length(fileNames):-1:1
    existingFiles(fileNameIndex) = exist(fileNames{fileNameIndex},'file') == 2;
end

if(~any(existingFiles))
    overwrite = 1;
    return;
end

overwrite = strcmp(questdlg('Tissue definition and/or computation results by this name already exist. Delete existing files?','Overwrite prompt','Yes','No, abort','Yes'),'Yes');

if ~overwrite
    fprintf('Aborted without saving data.\n');
    return;
end

for fileName = fileNames(existingFiles)
    delete(fileName{1});
end

return