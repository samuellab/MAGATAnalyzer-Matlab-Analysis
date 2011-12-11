function [dx, dy] = avggridspacing(x,y)
%average spacing between grid points, given points in grid, in no
%particular order
[y,I] = sort(y); x = x(I);
jmag = mean(diff(y)) + std(diff(y));
jump = [0 find(diff(y) > jmag)];
for j = 1:(length(jump) - 1)
    xx = x((jump(j) + 1):(jump(j+1)));
    yy = y((jump(j) + 1):(jump(j+1)));
    [xx,I] = sort(xx);
    yy = yy(I);
    row(j).xc = xx; %#ok<AGROW>
    row(j).yc = yy; %#ok<AGROW>
    dx(j) = median(diff(row(j).xc));
    yavg(j) = mean(row(j).yc);
end

dx = median(dx(isfinite(dx)));
dy = median(diff(yavg(isfinite(yavg))));