classdef model
    %MODEL Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        G MCmatlab.geometry
        MC MCmatlab.monteCarloSimulation
        FMC MCmatlab.fluorescenceMonteCarloSimulation
        HS MCmatlab.heatSimulation
    end
    
    methods
        function obj = model()
            %MODEL Construct an instance of this class
            %   Detailed explanation goes here
            
            obj.G = MCmatlab.geometry;
            obj.MC = MCmatlab.monteCarloSimulation;
            obj.FMC = MCmatlab.fluorescenceMonteCarloSimulation;
            obj.HS = MCmatlab.heatSimulation;
            
        end
        
        function obj = clearMCmatlabModel(obj, type)
            %clearMCmatlabModel Summary of this method goes here
            %   Detailed explanation goes here
            
            switch type
                case "G"
                    obj.G = MCmatlab.geometry;
                case "MC"
                    obj.MC = MCmatlab.monteCarloSimulation;
                case "FMC"
                    obj.FMC = MCmatlab.fluorescenceMonteCarloSimulation;
                case "HS"
                    obj.HS = MCmatlab.heatSimulation;
            end
        end
        
        obj = runMonteCarlo(obj, varargin)
       
        function plotMCmatlab(obj, varargin)
            %plotMCmatlab Summary of this method goes here
            %   Detailed explanation goes here
            
            if nargin == 1
                plotMCmatlabMC(obj)
                return
            end
            switch varargin
                case "fluorescence"
                    plotMCmatlabMC(obj, "fluorescence")
                case "geometry"
                    plotMCmatlabGeom(obj)
                case "heatSimulation"
                    plotMCmatlabHeat(obj)
            end
        end
        
        obj = defineGeometry(obj)
        
        obj = simulateHeatDistribution(obj)
        
        plotMCmatlabMC(obj, varargin)
        plotMCmatlabGeom(obj)
        plotMCmatlabHeat(obj)

    end
    
    methods (Access = private)
        finiteElementHeatPropagator(T,Omega,heatSimParameters)
        finiteElementHeatPropagator_CUDA(T,Omega,heatSimParameters)
        
        getMediaProperties_funcHandles(model,simType)
        getOpticalMediaProperties(model,simType)
        getThermalMediaProperties(model)
        
        plotMediaProperties(nFig,model,simType)
        plotVolumetric(nFig,xraw,yraw,zraw,Mraw,varargin)
        updateVolumetric(h_f,M)
        inferno(m)
        
        MCmatlab(model,simType)
        MCmatlab_CUDA(model,simType)
    end
    
end

