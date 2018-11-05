classdef Track < handle
    %Track of points loaded from an experiment (.bin) file
    %Tracks can be segmented into runs and reorientations
    %Also includes functions for getting derived information about track
    %(e.g. velocity, path curvature) and plotting tracks
    
    properties(Transient = true)
        expt; %experiment this track came fro
    end
    
  %  properties (AbortSet = true, Transient = false)
        
%   end
    
    properties
        pt; %array of track points
        dr = DerivationRules(); %rules for interpolating, smoothing, and differentiating 
        so = WormSegmentOptions(); %rules for segmenting tracks into runs
        npts = 0; %total number of points
        nt = 1; %number of tracks on disk that make up this track
        locInFile = 0; %the position (in bytes from beginning of file) of the track in the file
        startFrame = 0; %the index of the first frame in which this track appears
        endFrame = 0; %the index of the last frame in which this track appears
                isrun = []; %isrun(j) is true if eti(j) is during a run
        iscollision = []; %iscollision(j) is true if eti(j) is during a close pass with another track
        
        dq = []; %derived quantities 
        saveData = []; %store data from transient events to be saved here
    end
    properties (Transient = true, AbortSet = true) %no longer save segmentation
        run = []; %periods of forward motion between reorientations
        reorientation = []; %periods of direction change between runs
    end
    properties (Dependent = true)
        trackNum; %track.expt.track(trackNum) = track
    end
    
    methods %constructor
        function t = Track(varargin)
            %constructor
            if (nargin >= 1 && isa(varargin{1}, 'Track'))
                t.clone(varargin{1});
            end
        end
    end
    methods %access methods
        function set.reorientation(obj, value)
            if (~isempty(value) && isa(value, 'TrackPart'))
                [value.track] = deal(obj);
            end
            obj.reorientation = value;
        end
        function set.run(obj, value)
            if (~isempty(value) && isa(value, 'TrackPart'))

                [value.track] = deal(obj);
            end
            obj.run = value;
        end
        function val = get.trackNum(obj)
            if (isempty(obj) || isempty(obj.expt) || ~isa(obj.expt, 'Experiment') || isempty(obj.expt.track) || ~isa(obj.expt.track, 'Track'))
                val = [];
                return;
            end
            val = find(obj.expt.track == obj);
        end
        function set.trackNum(obj, value)
            return;
        end
    end
    
    methods(Static)
         track = fromFile(fid, ptType, loadImageByIndex, loadContour, camcalinfo, minpts);
         varargout = validDQName (varargin);
         track = fromJava(jTr, ptType, trInd, loadImageByIndex, loadContour, camcalinfo, minpts);
         function track = fromMatFile(fname)
            track = loadObjectTypeFromMatFile(fname, 'Track');
         end
    end
    
    methods(Static, Access = protected)
        varargout = nameInList (list, varargin);       
    end
    
    methods 
        plotPath(track, pathType, linetype, varargin);
 %       plotSegmentation(track, varargin);
        h = plotColorPath(track, fieldName, varargin);
        h = plotFields (track, xfield, yfields, varargin);
        addTime(track, indx, elapsedTime);
        calculateDerivedQuantity(track, quantityName, recalculate);
        recalculateDerivedQuantities(track, varargin); %recalculate all already derived quantities
        qvec = getDerivedQuantity(track, quantityName, recalculate, varargin);
        qvec = getSubFieldDQ (track, subfield, quantityName, varargin); 
        c = precedes(track1, track2, maxFrameDiff, maxDist)
        [pt, ind, dist] = nearestPoint(track, loc)
        inds = indsAtTime(track, elapsedTime);
 %       segmentTrack (track, wormSegmentOptions)
        merge(track, track2)
        track2 = split(track, splitnum);
        playMovie(track, varargin)
        clickMovie(track, varargin)
        movieOps = makeMovie(track, movieOps, varargin);
        
        qvec = fieldAtTime(track, quantityName, elapsedTime);
        addGlobalQuantity (track, quantityName, time, value);
        deleteMe = trim (track, timerange, validrect);
        [ps, f, peakfs] = powerSpectrum(track, quantityName, timeInterval, varargin);
        [xc, np, tx] = crosscorrelate (track, fieldname1, fieldname2, varargin);
        [ac, np, tx] = autocorrelate (track, fieldname, varargin);
        str = getReport(track, startInd, endInd, varargin);
        tf = passesThroughBox (track, rect, timeInterval, varargin);
        function toMatFile(track, fname)
            save(fixFileNameWin(fname), 'track');
        end
        
        
        
    end
    
    
    
    methods(Access = protected)
        function mo = makeMovieTrackSpecific(track, mo, pt, iind) %#ok<MANU>
            return; % put in specific stuff for maggots, worms, etc. in appropriate tracks
        end
        function clone(tr, oldtr)
            if (~isa(oldtr, 'Track'))
                return;
            end
            fl = intersect(fieldnames(tr),fieldnames(oldtr));
            for j = 1:length(fl)
                tr.(fl{j}) = oldtr.(fl{j});
            end
        end
    end
           
    
end

