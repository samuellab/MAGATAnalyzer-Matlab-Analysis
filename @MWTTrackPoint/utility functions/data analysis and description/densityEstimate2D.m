function [density, xx, yy] = densityEstimate2D (data, xrange, yrange, bwx, bwy, M)
%function [density, xx, yy] = densityEstimate2D (data, xrange, yrange, bwx, bwy, M)
%
%data is N x 2
if (size(data,2) > size(data,1))
    warning('transposing data to make N x d');
    data = data';
end
dx = abs(diff(xrange))/M;
dy = abs(diff(yrange))/M;

padx = ceil(6*bwx/dx)+1;
pady = ceil(6*bwy/dy)+1;
pad = max(padx,pady);
MM = M + 2*pad;
data(:,1) = (data(:,1)-min(xrange-dx*pad))/(abs(diff(xrange)) + 2*dx*pad);
data(:,2) = (data(:,2)-min(yrange-dy*pad))/(abs(diff(yrange)) + 2*dy*pad);

h = ndhist(data, MM);
hs = conv2(gaussKernel(bwy/dy), gaussKernel(bwx/dx), h, 'same');
density = hs(pad + (1:M), pad + (1:M));
xaxis = min(xrange) - 0.5 * dx + (1:M)*dx;
yaxis = min(yrange) - 0.5 * dy + (1:M)*dy;
[xx,yy] = meshgrid(xaxis,yaxis);

