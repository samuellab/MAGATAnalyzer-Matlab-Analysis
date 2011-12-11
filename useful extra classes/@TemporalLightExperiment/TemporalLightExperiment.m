classdef TemporalLightExperiment < Experiment
    %an experiment that has the lights turned on and off; a la carter
    
    properties
        cycleTime = []; %frames per half period
    end
    
    methods
        cycleTime = detectCycleTime (expt, varargin);
        addCycleTime (expt, varargin);
    end
    
    methods %constructor
        function expt = TemporalLightExperiment(varargin)
            expt = expt@Experiment(varargin{:});
            if (nargin > 1)
                expt.cycleTime = varargin{2}; 
            else
                expt.detectCycleTime();
            end
            if (nargin > 0)
                expt.addCycleTime();
            end
        end
    end
end

