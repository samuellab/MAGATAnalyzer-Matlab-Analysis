classdef ImTrackPoint < TrackPoint
    %Augments TrackPoint by including an excerpted image
    
    
    properties
        imOffset = int16([0 0]); %location of imagedata origin from camera origin (always in pixels)
        imData = uint8([]); %matrix representing image
    end
    methods
        tp = fromFile (tp, fid, loadIm, loadContour, camcalinfo)
        drawTrackImage(tp, camcalinfo, varargin)
    end
    
    methods %constructor
        function tp = ImTrackPoint(varargin)
            tp = tp@TrackPoint(varargin{:});            
            if (nargin >= 1) && (isa(varargin{1}, 'ImTrackPoint'))
                op = varargin{1};
                flist = fieldnames(tp);
                for j = 1:length(flist)
                     tp.(flist{j}) = op.(flist{j});
                end
            end
        end%trackpoint
    end
    
end

