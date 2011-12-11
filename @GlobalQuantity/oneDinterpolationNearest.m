function yout = oneDinterpolationNearest(xin, xData, yData)
% does a 1-D lookup table with nearest neighbor interpolation
% function yout = oneDinterpolation(xin, xData, yData)

%reshape y to have correct dimension
sz = size(yData);
n = length(xData);

ind = find(sz == n);
ind2 = find(sz ~= n);

if (~any(ind))
    errmsg = ['yData and xData have incommensurate sizes: size(yData) = ' mat2str(sz) ...
        ' and size xData = ' mat2str(size(xData))];
    error('GERSHOW:GQ01', errmsg);
end

order = [ind ind2];
yData = permute(yData, order);
size(yData);

%make sure xin is a column vector
if (size(xin,2) == 1)
    xin = xin';
end
s = warning('off','all');
%interpolate
yout = interp1 (double(xData), double(yData), xin, 'nearest', NaN);
warning(s);
end