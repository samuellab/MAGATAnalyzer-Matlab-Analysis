classdef WormReorientation < TrackPart
    %a worm reorientation, consisting of multiple sharp turns
    
    properties
        sharpTurn;
        numTurns;
        thetaIn;
        thetaOut;
        dTheta;
        turnsequence;
    end
    properties (Dependent = true)
        prevRun;
        nextRun;
        prevDir;
        nextDir;
    end
    
    methods
        function nr = get.nextRun(obj)
            nr = obj.getAdjacent('next', 'run');
        end
        function pr = get.prevRun(obj)
            pr = obj.getAdjacent('prev', 'run');
        end
        function nd = get.nextDir(obj)
            nd = obj.thetaOut;
        end
        function pd = get.prevDir(obj)
            pd = obj.thetaIn;
        end
    end
    methods
        addTurn(reo, st);
        calculateMetrics(reo);
        st = firstTurn(reo); 
        str = getReport(reo, varargin);
        tf = turnsequenceEquals(reo, turnsequence);
    end
    
end

