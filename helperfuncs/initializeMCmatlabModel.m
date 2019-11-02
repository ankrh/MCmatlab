function model = initializeMCmatlabModel()

model = struct;

model = clearMCmatlabModel(model,'G');
model = clearMCmatlabModel(model,'MC');
model = clearMCmatlabModel(model,'FMC');
model = clearMCmatlabModel(model,'HS');
end