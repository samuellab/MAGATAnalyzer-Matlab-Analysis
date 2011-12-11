function [area, areas, threshold] = optimumArea(track, varargin)
%function area = optimumArea(track, varargin)
%
%optional args
%step = length(track.pt)/200;
%minThresh = 1;
%maxThresh = 255;
%minArea = 5;
%maxArea = 500;
%scale = -1;
%blur = true;
%showGraphs = false;

blur = false;
scale = -1;
step = [];
varargin = assignApplicable(varargin);
track.expt.openDataFile;

pt = [track.pt];
if (isempty(step))
    step = ceil(length(pt)/250);
end
if ~(isa(pt, 'ImTrackPoint'))
    error('Need a track of image points');
end
inds = 1:step:length(pt);
areas = zeros(size(inds));
threshold = zeros(size(inds));

for j = 1:length(inds)
    if isempty(pt(inds(j)).imData)
        pt2 = track.expt.reloadPoint(pt(inds(j)));
        im = pt2.imData;
    else
        im = pt(inds(j)).imData;
    end
    if (scale > 0)
        im = imresize(im, scale);
    end
    if (blur)
        imnb = im;
        scale = max(scale, 1);
        im = blurim(im,scale);
    end
    [t,a] = optimumThreshold(im, varargin{:});
    if (blur)
        a = sum(imnb(:) > t);
    end
    areas(j) = a/scale^2;
    threshold(j) = t;
end
area = median(areas);