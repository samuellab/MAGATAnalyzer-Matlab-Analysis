function [xout, yout, inds] = nearestNPoints (x0, y0, x, y, npts)
%function [xout, yout] = nearestNPoints (x0, y0, x, y, npts)

d = (x-x0).^2 + (y - y0).^2;
[~,I] = sort(d);
xout = x(I(1:npts));
yout = y(I(1:npts));
inds = I(1:npts);