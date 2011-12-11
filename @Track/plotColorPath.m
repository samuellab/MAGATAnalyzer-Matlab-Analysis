function h = plotColorPath(track, fieldName, varargin)
% plots track path with intensity given by track.dq.fieldName
% function plotColorPath(track, fieldName, varargin)
% plots track path with intensity given by track.dq.fieldName
%
% outputs:
%   H: (optional) handles to all points in the line
% inputs:
%   TRACK < Track
%   FIELDNAME: the name of the field to plot
%   VARARGIN: parameter/value pairs
%   defaults that can be overwritten:
%   pathType = 'sloc' {'loc', 'iloc'}
%   zdata = track.dq.(fieldName) or empty if fieldName is empty
%       to put in different zdata, leave fieldName empty and add pair
%       'zdata', zvalue
%   inds = all points
%       if zdata has the same size as path, we plot path(:,inds) vs.
%       zdata(inds)
%       else if inds has the same size as zdata diff from path, we plot 
%       path(:,inds) vs. zdata
%   zrange = range over which the zdata is scaled (default min(z) max(z))
%   cmap = colormap to use (default jet(256))


pathType = 'sloc';
if (~isempty(fieldName))
    zdata = track.getDerivedQuantity(fieldName);
else
    zdata = [];
end
inds = [];
zrange = [];
cmap = jet();
varargin = assignApplicable(varargin);

path = track.getDerivedQuantity(pathType);

if (isempty(inds))
    inds = 1:length(path);
end
if (length(zdata) == length(path))
    zdata = zdata(inds);
else
    if (length(zdata) ~= length(inds))
        disp('length of zdata doesn''t match either length of path or length of inds');
        disp(['length(inds) = ' num2str(length(inds)) ' length(zdata) = ' num2str(length(zdata)) ' length(path) = ' num2str(length(path))]);
        return;
    end
end
    
if (nargout > 0)
    h = plotColorLine (path(1,inds), path(2,inds), zdata, cmap, zrange, varargin{:});
else
    plotColorLine (path(1,inds), path(2,inds), zdata, cmap, zrange, varargin{:});
end
    