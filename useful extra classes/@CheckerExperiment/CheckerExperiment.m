classdef CheckerExperiment < Experiment
    %A class that includes some special functions for setting up the
    %analysis of one of Ashley and Liz's famous checkerboard experiments
    %also, cures cancer
    
    properties
        structvar; %There is no god but structvar and structvar is his prophet
        boundaryThickness = 30;
        imsize = [1944 2592];
    end
    
    methods
        addSpatialInfo(expt, varargin);
        addBoundary(expt, dirim,distim,lightim);
        setBoundaryThickness(expt, thickness);
    end
    
    methods %constructor
        function expt = CheckerExperiment(varargin)
            expt = expt@Experiment(varargin{:});
            if (length(varargin) > 1 && ischar(varargin{2}))
                try
                    d = dir(varargin{2});
                catch %#ok<CTCH>
                    return;
                end
                if (length(d) == 1)
                    load(varargin{2}, 'structvar');
                    expt.structvar = structvar; %#ok<CPROP,*PROP>
                end
            end
            
        end
    end
    methods (Static)
        [dirim,distim,interior] = generateBoundaryImage(imsize, boundaryThickness, structvar);
    end
end

