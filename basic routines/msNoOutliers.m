function [m, s] = msNoOutliers (y)
%function [m, s] = msNoOutliers (y)

m = mean(y);
s = std(y);
j = 0;
dm = m;
ds = s;
while (j < 100 && (abs(dm/m) > 0.001) && (abs(ds/s) > 0.001))
    inds = find(abs(y - m) < 1.5*s);
    oldm = m;
    olds = s;
    m = mean(y(inds));
    s = std(y(inds));
    dm = m - oldm;
    ds = s - olds;
    j = j+1;
end

