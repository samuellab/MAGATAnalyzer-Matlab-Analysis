classdef Run < TrackPart
    %UNTITLED7 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        startTheta;
        endTheta;
        meanTheta;
        euclidLength;
        pathLength;
        runTime;    
    end
    properties (Dependent)
        nextRun;
        previousRun;    
        prevReorientation;
        nextReorientation;
    end
    methods
        function nr = get.nextRun(obj)
            nr = obj.getAdjacent('next', 'run');
        end
        function nr = get.nextReorientation(obj)
            nr = obj.getAdjacent('next', 'reorientation');
        end
        function pr = get.previousRun(obj)
            pr = obj.getAdjacent('prev', 'run');
        end
        function pr = get.prevReorientation(obj)
            pr = obj.getAdjacent('prev', 'reorientation');
        end
    end    
    
    %{
    methods
        function set.previousRun (obj, value)
            if (~isempty(value) && isa (value, 'Run'))
                value.nextRun = obj;
            end
            obj.previousRun = value;
        end
        function set.prevReorientation (obj, value)
            if (~isempty(value) && isfield (value, 'nextRun'))
                value.nextRun = obj;
            end
            obj.prevReorientation = value;
        end
        function set.nextReorientation (obj, value)
            if (~isempty(value) && isfield (value, 'prevRun'))
                value.nextRun = obj;
            end
            obj.nextReorientation = value;
        end
    end
    %}
    methods
        calculateMetrics(run);
        draw(run, varargin);
        str = getReport(run, varargin);
    end
    methods %constructor
        function run = Run(varargin)
            %run = Run(track, startInd, endInd)
%            dbstack
            arglist = {'track', 'startInd', 'endInd'};
            if (nargin > 0)
                for j = 1:min(nargin, length(arglist))
                    run.(arglist{j}) = varargin{j};
                end
                if (nargin >= length(arglist))
                    %got everything we wanted, so calculate away
                    run.calculateMetrics();
                end
            end
        end
    end
end

