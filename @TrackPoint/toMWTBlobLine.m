function str = toMWTBlobLine(tp, camcalinfo, varargin)
% function str = toMWTBlobLine(tp, camcalinfo, varargin)
% produces a line of text corresponding to a line in the mwt blob file
%
% tp < MaggotTrackPoint
% str - character array

if (~exist('camcalinfo', 'var')), camcalinfo = []; end
varargin = assignApplicable(varargin);

data = zeros(10, 1);

data(1) = round(tp.ind);
data(2) = tp.et;
if (~isempty(camcalinfo))
    data(3:4) = camcalinfo.camPtsFromRealPts(tp.loc);
else
    data(3:4) = tp.loc;
end
data(5) = round(tp.area);
c = [tp.cov(1) tp.cov(2); tp.cov(2) tp.cov(3)];

[V,D] = eig(c);
[~,I] = sort(sum(D));
data(6:7) = V(:,I(2))*D(I(2),I(2));
data(8) = D(I(1),I(1));

%MWT does a fit of the maggot to a line, then gives the maximum length of a line
%parallel to that fit line with endpoints on the contour, and a similar length for the perpendicular line
%in order to avoid loading the pixel values from disk again, we will just
%use the longest line within the contour and the longest line
%perpendicular to that one

if (~isempty(camcalinfo))
    cpts = round(camcalinfo.camPtsFromRealPts(tp.contour));
else
    cpts = tp.contour;
end
x = repmat(cpts(1,:), size(cpts,2), 1);
y = repmat(cpts(2,:), size(cpts,2), 1);
dd = (x-x').^2 + (y-y').^2; 
[dm,I] = max(dd(:));
%longest distance
data(9) = sqrt(dm);
[I,J] = ind2sub(size(dd), I);

%perpendicular vector
vp = [-diff(cpts(2, [I J])) diff(cpts(1,[I J]))];
vp = vp./sqrt(sum(vp.^2));
dd = abs(vp(1)*(x-x') + vp(2)*(y-y'));
data(10) = max(dd(:));
str = sprintf('%8g\t', data);
