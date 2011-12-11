classdef MaggotReorientation < TrackPart
    % periods of reorientation between runs; in our publications, a
    % reorientation is called a "turn"
    % reorientations have zero or more headSwings ("head sweeps" in
    %   publications)
    
    properties
        numHS = 0; %number of headswings in reorientation
        headSwing; %pointer to headswings in reorientation
        prevDir = NaN; %heading direction at end of previous run 
        nextDir = NaN; %heading direction at start of next run
    end
    properties(Dependent = true)
        prevRun;
        nextRun;
    end
    
    
    methods %constructor
        function reo = MaggotReorientation (track, headswings, prevRun, nextRun)
            %if no headswings (headswings = []), then pass prevRun, nextRun
            %to specify interval
            switch nargin
                case 0
                    return;
                case 1
                    reo.track = track;
                    return
                otherwise
                    reo.track = track;
                    if (~isempty(headswings))
                        for j = 1:length(headswings)
                            headswings(j).num = j;
                        end
                    end
                    reo.headSwing = headswings;
                        
                    existsAndDefault ('prevRun', []);
                    existsAndDefault ('nextRun', []);
                    reo.calculateMetrics (prevRun, nextRun);
            end
        end
    end
    
    methods
        function nr = get.nextRun(obj)
            nr = obj.getAdjacent('next', 'run');
        end
        function pr = get.prevRun(obj)
            pr = obj.getAdjacent('prev', 'run');
        end
    end
    
    
    methods 
        calculateMetrics(reo, prevRun, nextRun);
        draw(reo, varargin)
    end
end

