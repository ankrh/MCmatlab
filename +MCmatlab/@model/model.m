classdef model
    %MODEL Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        G MCmatlab.geometry
        MC MCmatlab.monteCarloSimulation
        FMC MCmatlab.monteCarloSimulation
        HS MCmatlab.heatSimulation
    end
    
    methods
        function obj = model()
            %MODEL Construct an instance of this class
            %   Detailed explanation goes here
            
            obj.G = MCmatlab.geometry;
            obj.MC = MCmatlab.monteCarloSimulation;
            obj.FMC = MCmatlab.monteCarloSimulation;
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
                    obj.FMC = MCmatlab.monteCarloSimulation;
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
end

