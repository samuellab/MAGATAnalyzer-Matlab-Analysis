function drawContourAndHead(pt, varargin)
%function drawContourAndHead (pt, varargin)
%@MaggotTrackPoint
%
%optional arguments and defaults
%
%offset = [0;0];
%contourColor = 'k-';
%headMarker = 'g*';
%tailMarker = 'rh';


offset = [0;0];
contourColor = 'k-';
spineColor = 'y-';
headMarker = 'g*';
tailMarker = 'rh';
varargin = assignApplicable(varargin);
h = pt.head + offset;
m = pt.mid + offset;
t = pt.tail + offset;
c = pt.contour + repmat(offset, 1, length(pt.contour));
c(:,end+1) = c(:,1); 


plot (h(1),h(2),headMarker, t(1),t(2),tailMarker, [t(1),m(1),h(1)], [t(2),m(2),h(2)], spineColor, c(1,:), c(2,:), contourColor, varargin{:});
