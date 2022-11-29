function plotAzFz(model)
  MCmatlab.NdimSliderPlot(model.G.dx*model.G.dy*sum(sum(model.MC.NA ,1),2),'nFig',41,'axisLabels',{'','','z [cm]','Az [cm^{-1}]'},'axisValues',{[],[],model.G.z});
  MCmatlab.NdimSliderPlot(model.G.dx*model.G.dy*sum(sum(model.MC.NFR,1),2),'nFig',42,'axisLabels',{'','','z [cm]','Fz [-]'},'axisValues',{[],[],model.G.z});
end