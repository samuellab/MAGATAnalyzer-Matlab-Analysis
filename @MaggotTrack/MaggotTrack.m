classdef MaggotTrack < Track
    %MaggotTrack extends Track 
    %MaggotTrack adds the headSwing field (keeps track of when maggot is
    %   turning its head to the side)
    %Several Track methods (e.g. segmentTrack) are overwritten or extended
    %   here
    %
    %Typical end users will not instantiate their own MaggotTracks but
    %   will find them in experiments they have loaded from disk
    
    properties (Transient = true, AbortSet = true)
        headSwing;
    end
    
    methods %constructor
        function mt = MaggotTrack (varargin)
            mt = mt@Track(varargin);
            mt.dr.smoothTime = 0.5; %maggots move faster & need less smoothing than worms
            mt.so = MaggotSegmentOptions();
            if ((nargin >= 1) && isa(varargin{1}, 'MaggotTrack'))
                mt.clone(varargin{1});
            end
        end
    end %constructor
    
    methods %access
        function set.headSwing(obj, value)
            if (~isempty(value) && isa(value, 'TrackPart'))
                [value.track] = deal(obj);
            end
            obj.headSwing = value;
        end
    end
    
    methods(Static)
        varargout = validDQName (varargin)
    end
    methods 
        recalculateDerivedQuantities(track, varargin);   %recalculate already derived quantities
        calculateDerivedQuantity(track, quantityName, recalculate); 
        %setSegmentSpeedsByPercentile (track, maggotSegmentOptions, stopPctl, startPctl); %nix 
        segmentTrack (track, maggotSegmentOptions); %segments the track into runs and turns (with head sweeps)
        setSegmentSpeeds (track, mso);
        plotSegmentation (track, varargin); %plots the segmented tracks run and turns with annotation
        fixHTOrientation(track, varargin); 
        prettyMovie(track, varargin); %play a pretty movie
        
        toMWTTrackFile (track, filestub, varargin); %writes a .blob file for 1 track
        str = toMWTTrackString (track, varargin); %writes the appropriate string for a .blobs file
        shapeModel = getTypicalMaggotShape(track, varargin);
        
        playMovie_BehaviorTriggered(track, fieldName,fieldDescription, displacementAxis, triggeredSum, triggeredInd, figTitle, varargin);
        loadJavaMasks(track, jTr);
        
    end
    
        
    
    
end

