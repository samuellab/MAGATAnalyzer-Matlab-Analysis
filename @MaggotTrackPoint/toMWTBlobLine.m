function str = toMWTBlobLine(tp, camcalinfo, varargin)
% function str = toMWTBlobLine(tp, camcalinfo, varargin)
% produces a line of text corresponding to a line in the mwt blob file
%
% tp < MaggotTrackPoint
% str - character array

if (nargin < 2), camcalinfo = []; end
varargin = assignApplicable(varargin);

data = zeros(10, 1);

data(1) = round(tp.ind) + 1; %MWT starts with frame #1, we start with frame 0
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
str1 = [sprintf('%d ', data(1)) sprintf('%8g ', data(2:end))];

if (isempty(tp.spine))
    tp.spine = repmat(tp.loc, 11);
end
if ~(isempty(camcalinfo))
    sp = camcalinfo.camPtsFromRealPts(tp.spine);
else
    sp = tp.spine;
end

sp = interp1(find(all(isfinite(sp))), sp(:,all(isfinite(sp)))', linspace(length(sp),1, 11))'; %reverse spine direction so head is first (MWT spinesForward plugin result)

sp(1,:) = round(sp(1,:) - data(3));
sp(2,:) = round(sp(2,:) - data(4));
str2 = sprintf('%d ', sp(:));

%create chain code for contour
cpts = removeSmallLoops(cpts);
dx = round(diff(cpts(1,[1:end 1])));
dy = round(diff(cpts(2,[1:end 1])));

while (any(abs(dx) > 1)) 
    I = find(abs(dx) > 1, 1, 'first');
    n = abs(dx(I));
    dx = [dx(1:(I-1)) repmat(sign(dx(I)), 1, n) dx((I+1):end)];
    dy = [dy(1:I) zeros(1,n-1) dy((I+1):end)];
end

while (any(abs(dy) > 1)) 
    I = find(abs(dy) > 1, 1, 'first');
    n = abs(dy(I));
    dy = [dy(1:(I-1)) repmat(sign(dy(I)), 1, n) dy((I+1):end)];
    dx = [dx(1:I) zeros(1,n-1) dx((I+1):end)];
end
[dx,dy] = convert8Connectedto4Connected(dx,dy);
if (length(dx) < 10)
    disp('bad contour!'); 
    tp.ind
end
str3 = sprintf('%d %d %d %s', round(cpts(1,1)), round(cpts(2,1)), length(dx), toKerrChainCode(dx,dy));

str = [str1 ' % ' str2 ' %% ' str3];