classdef HeadSwing < TrackPart
    %what happens when a Maggot swings its head
    %it's a party, that's what!
    
    properties
        maxInd; %index where maxTheta is found; i.e. the peak of the headsweep
        maxTheta = 0;
        headDir = 0; 
        prevDir = NaN;
        nextDir = NaN;
        tailDir;
        sign;
        accepted;
        num = 0; %which head sweep of the reorientation
    end
    properties (Dependent = true)
        prevRun;
        nextRun;
    end
    methods
        function nr = get.nextRun(obj)
            nr = obj.getAdjacent('next', 'run');
        end
        function pr = get.prevRun(obj)
            pr = obj.getAdjacent('prev', 'run');
        end
    end
    
    methods %constructor
        function hs = HeadSwing(varargin)
            %run = Run(track, startInd, endInd)
            arglist = {'track', 'startInd', 'endInd'};
            if (nargin > 0)
                for j = 1:min(nargin, length(arglist))
                    hs.(arglist{j}) = varargin{j};
                end
                if (nargin >= length(arglist))
                    %got everything we wanted, so calculate away
                    hs.calculateMetrics();
                end
            end
        end
    end %constructor
    methods
        calculateMetrics(headSwing);
        draw(headSwing, varargin)
        illustrateHeadSwing (hs, varargin);
    end
end

