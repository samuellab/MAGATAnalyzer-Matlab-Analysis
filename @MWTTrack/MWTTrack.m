classdef MWTTrack < MaggotTrack
    %Extension of MaggotTrack for Multi Worm Tracker blob files
    %main difference is how to load a file
    
    properties
        fname;
    end
    methods
        fixHTOrientation(track, varargin);
        markHTInvalid(track, thresh, varargin);
        valid = removeCollisionPoints(track, maxArea, varargin);        
        calculateDerivedQuantity(track, quantityNames, recalculate);
    end
    methods %constructor
        function mt = MWTTrack (varargin)
            mt = mt@MaggotTrack(varargin);
            mt.dr.interpTime = 0.033; %video taken at ~30 fps
            mt.dr.smoothTime = 0.25; %try less smoothing & see
            mt.dr.derivTime = 0.25; %1 pt derivative
            mt.so = MaggotSegmentOptions();
            mt.so.speed_field = 'smoothSpeed';
            if ((nargin >= 1) && isa(varargin{1}, 'MaggotTrack'))
                mt.clone(varargin{1});
            end
        end
    end %constructor
    methods (Static)
        mt = fromFile(fname, camcalinfo); %or fid
        varargout = validDQName (varargin);
    end
    
end

