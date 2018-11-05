function [m, s] = msNoOutliers (y, sigscale)
%function [m, s] = msNoOutliers (y, sigscale)

m = mean(y);
s = std(y);
j = 0;
dm = m;
ds = s;
if (nargin < 2)
    sigscale = 1.5;
end
while (j < 100 && (abs(dm/m) > 0.001) && (abs(ds/s) > 0.001))
    inds = find(abs(y - m) < sigscale*s);
    oldm = m;
    olds = s;
    m = mean(y(inds));
    s = std(y(inds));
    dm = m - oldm;
    ds = s - olds;
    j = j+1;
end

