classdef TrackPart < handle
    %Basic superclass for anything that labels a section of track
    %for instance, runs, reorientations, headswings, potty breaks
    
    properties(Transient = true, AbortSet = true)
         track;
    end
    
    properties
       
        startInd = 0;
        endInd = 0;
        inds = [];
        valid = true; %note, not the same as isvalid
    end
    
    methods
        qvec = getDerivedQuantity (tp, field, varargin);
        [qv, datamatrix] = averageDerivedQuantity (tp, field, centerpos, offsetinds, varargin);
        
        playMovie (tp, varargin);
        movieOps = makeMovie(track, movieOps, varargin);
        prettyMovie (tp, varargin);
        tf = containsIndex(tp, ind); 
        str = getReport(tp, varargin);
        h = plotFields (tp, xfield, yfields, varargin);
        plotPath(tp, pathType, linetype, varargin);
    end
    methods (Access = protected)
        myMovie (tp, mtype, varargin);
        
        function atp = getAdjacent(tp, direction, fieldname)
            if (isempty(tp.track) || ~isa(tp.track, 'Track'))
                disp ('trackpart doesn''t point to track');
                atp = [];
                return;
            end
            t = tp.track;
            f = [t.(fieldname)];
            if (~isa(f, 'TrackPart'))
                %disp ('desired field isn''t a trackpart');
                atp = [];
                return
            end
            [~,I] = sort([f.startInd]);
            f = f(I);
            if (direction(1) == 'n')
                ind = find([f.endInd] > tp.endInd, 1, 'first');
                atp = f(ind);
            else
                ind = find([f.startInd] < tp.startInd, 1, 'last');
                atp = f(ind);
            end
                
        end
            
                
    end
    
end

