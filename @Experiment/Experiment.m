%Experiment
%A set of tracks extracted from the same experiment
%Each experiment is loaded from a single .bin file created by the
%track extraction software


classdef Experiment < handle
    %Experiment
    %A set of tracks extracted from the same experiment
    %Each experiment is loaded from a single .bin file created by the
    %track extraction software

    properties
        fname = ''; %the name of the .bin file representing the experiment
        timfname = ''; %the name of the timing file that tells when each frame happened
        camcalinfo = []; %camera calibration information; not tested at the moment
        elapsedTime = []; %the elapsed time from the start of the experiment for each frame
        %temperature = []; %??? nix
        
        dr = DerivationRules(); %rules for interpolating, smoothing, and differentiating
        so = WormSegmentOptions(); %rules for segmenting tracks into runs & reorientations
        globalQuantity; %used to add extra information (e.g. light, temperature) to the tracks
    end
    properties (AbortSet = true)
        track; %the set of tracks contained in the experiment; see Track, MaggotTrack
    end
    
    properties (Transient = true)
        fid = 0;
    end
    
    methods
        addtime (expt, timfname) %adds timing information to map frame number to seconds; should be done automatically on loading
        addGlobalQuantity(expt,varargin); %adds a global quantity (maps derived quantity, e.g. elapsed time) to another quantity (e.g. odor level)
        assignGlobalQuantities(expt,varargin); %assigns all global quantities to tracks; automatically called when quantity is added
        addStandardizedField(expt, fieldname, varargin); %standardized_field = (field - mean(field))/std(field) on a track by track basis
        addMetaDataFields(expt, timfname, varargin); %adds global quantities from meta data file (.mdat) 
        stitchTracks (expt, maxFrameInterval, maxDist, varargin); %stitches tracks together if they come from the same animal
        cleanTracks (expt, minFrames, minDist); %deprecated - use ESetCleaner
        qvec = gatherField(expt, fieldname, varargin); %gets all values of a fieldname
        qvec = gatherFieldInTracks(expt, fieldname, tracknums, varargin)%gets all values of 'fieldname' for a subset of tracks in expt
        qvec = gatherFromSubField(expt, subfield, fieldname, varargin); %gets all values of 'fieldname' from track parts (e.g. runs, reorientations)
        qvec = gatherSubField (expt, field, subfield, varargin); %gathers from a subfield of track; all track.(field).(subfield)
        [pt, track, trackind, ptind] = findNearestPoint (expt, x, y); %finds nearest point in the experiment to given point
        pt = reloadPoint (expt, pt); %reloads a point from disk, including image, if any
        track2 = reloadTrack(expt, track); %reloads a track from disk, including images, if any
        
        openDataFile(expt); %opens the data file the experiment was originally loaded from, for reading
        closeDataFile(expt); %closes the data file, if open
        calculateDerivedQuantity(expt, quantityNames, recalculate); %assigns experiments derivation rules to all tracks, then calculates derived quantity
        
        segmentTracks(expt, segmentOptions); %segments all tracks in experiment; rarely used (eset.executeTrackFunction('segmentTrack') is more typical call)
        varargout = executeTrackFunction(expt, func, varargin); %executes the function func for all tracks in expt
        result = evaluateTrackExpression(expt, expression); %evaluates the expression for all tracks in expt
        detectCollisions(expt, maxdist); %marks the iscollision field of each track if it comes too close to another track
        [inds,why] = detectPossibleSegmentationProblems(expt); %detects possible segmentation problems with maggot tracks (rarely used)
        im = makeWholeFrameImage (expt, ind, varargin); %creates a mosaic image of all the individual point images for a frame
        
        %trimTracks cuts out any part of the track outside
        %min(timerange),max(timerange) and removes any part of the track
        %from the point the track leaves validrect until the end of the
        %track; leave timerange or validrect empty to disable
        trimTracks(expt, timerange, validrect); %removes parts of tracks
        
        %pruneTracks removes completely any track that starts outside
        %starttimerange in elapsedTime, as well as any track that starts
        %outside startrect; leave startttimerange or startrect empty to disable
        pruneTracks(expt, starttimerange, startrect); %removes entire tracks
        
        [ps, f] = powerSpectrum(expt, quantityName, timeInterval, varargin); %produces the power spectrum for a qiven measured/derived quantity
        [xc, np, tx, nt] = crosscorrelate (expt, fieldname1, fieldname2, varargin);%cross correlates fieldname1 and fieldname2
        [ac, np, tx, nt] = autocorrelate (expt, fieldname, varargin);%autocorrelates fieldname
        updateTrackSegmentationOptions(expt); % changes all tracks segmentation options to be the same as experiment's
        setTimingByFrameDifference(expt, deltaT, overrideExisting); %sets the elapsedTime to be 0:1:numFrames * deltaT; another butt-saver that shouldn't need to be used
        
        addTonToff(expt, fieldname, ramptype, varargin); %adds _ton and _toff global quantities for a given temporal global quantity
    end
    
    methods %set methods
        function set.track(obj, value)
            if (~isempty(value) && isa (value, 'Track'))
                [value.expt] = deal(obj);
            end
            obj.track = value;
        end
        
    end
    
    methods %constructor
        function expt = Experiment(varargin)
            if (~isempty(varargin) && isa(varargin{1}, 'Experiment'))
                fn = intersect(fieldnames(varargin{1}), fieldnames(expt));
                for j = 1:length(fn)
                    expt.(fn{j}) = varargin{1}.(fn{j});
                end
            end
        end
    end
    
    
    methods(Static)
        expt = fromFile (fname, timfname, loadContour, camcalinfo, minTrackLength)
    end
end

