classdef MaggotTrackPoint < ImTrackPoint
    % MaggotTrackPoint extends ImTrackPoint by adding 
    % information about the head, tail, midline (spine), and contour
    
    properties
        targetArea = -1;
        threshold = -1;
        htValid = false;
        head = [NaN NaN];
        mid = [NaN NaN];
        tail = [NaN NaN];
        contour = [NaN NaN];
        spine = [NaN NaN];
    end
    
    methods
        tp = fromFile (tp, fid, loadIm, loadContour, camcalinfo);
        drawTrackImage(tp, camcalinfo, varargin);
        drawContourAndHead(pt, varargin);
        str = toMWTBlobLine(tp, camcalinfo, varargin);
        mtp2 = unCameraCalibrate (mtp, camcalinfo);
        tp = fromJava(tp, jTP, loadIm, loadContour, camcalinfo);
    end
    
    methods
        function tp = MaggotTrackPoint(varargin)
            tp = tp@ImTrackPoint(varargin{:});
            if (nargin >= 1) && (isa(varargin{1}, 'MaggotTrackPoint'))
                op = varargin{1};
                flist = fieldnames(tp);
                for j = 1:length(flist)
                    tp.(flist{j}) = op.(flist{j});
                end
            end
        end%MaggotTrackPoint
   
    end
    %}
end

