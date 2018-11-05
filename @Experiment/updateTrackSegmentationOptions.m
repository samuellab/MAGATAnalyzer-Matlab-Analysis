function updateTrackSegmentationOptions(expt, varargin)
% changes all tracks segmentation options to be the same as experiment's
% function updateTrackSegmentationOptions(expt, varargin)
%
% EXPT < experiment
if (length(expt) > 1)
    for j = 1:length(expt)
        expt(j).updateTrackSegmentationOptions(varargin{:});
    end
    return;
end
if (length(varargin) >= 1)
    expt.so = varargin{1};
end
[expt.track.so] = deal(expt.so);