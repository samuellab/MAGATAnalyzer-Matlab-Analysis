classdef WormTrack < Track
    %Tracks from c. elegans
    
    properties (Transient = true, AbortSet = true) %no longer save segmentation
        sharpTurn = []; %information about sharp turns in the track (worm only)
    end
    methods %constructor
        function t = WormTrack(varargin)
            %constructor
            if (nargin >= 1 && isa(varargin{1}, 'Track'))
                t.clone(varargin{1});
            end
        end
    end
    % save and load methods
    methods
        function track = saveobj(track)
            track = track.fillOutSaveData;
        end
        
    end
    methods (Static)
        function track = loadobj(track)
           % disp ('loadobj called on worm track');
           try
                if (~isempty(track.saveData))
                    track = track.restoreSaveData;
                end
           catch me
               disp(me.getReport);
           end
        end
    end
    
    methods %access methods
        function set.sharpTurn(obj, value)
            if (~isempty(value) && isa(value, 'TrackPart'))
                [value.track] = deal(obj);
            end
            obj.sharpTurn = value;
        end
    end
    methods 
        plotSegmentation(track, varargin);       
        segmentTrack (track, wormSegmentOptions, varargin);
        flagOmegaTurnsAndReversals (track, wormSegmentOptions, varargin);
        flagReorientations (track, wormSegmentOptions, varargin)

        playMovie(track, varargin)
        track = fillOutSaveData(track);
        track = restoreSaveData(track);
    end
    methods(Access = protected)
        mo = makeMovieTrackSpecific(track, mo, pt, iind) 
    end
    
end

