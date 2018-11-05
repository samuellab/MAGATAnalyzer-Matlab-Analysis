function yout = tri2Dinterpolation(xin, xData, yData)
% interpolates a function defined on a scattered set of 2D points
% function yout = tri2Dinterpolation(xin, xData, yData)
%
% uses delauney triangulation (TriScatteredInterp) to interpolate from
% scattered lists of points and values
%
% xin is a 2xM list of points 
% xData is a 2xN list of points
% yData is a kxN list of values 
% yData(:,j) = underlyingFunction(xdata(:,j));
% output is a kxM list of values
% output = underlyingFunction(xin)


%preallocate yout, preserving type
yout = repmat(yData(1), [size(yData,1),size(xin,2)]);
islog = islogical(yData);
xin = double(xin);
xData = double(xData);
yData = double(yData);


for k = 1:size(yData,1)
%     size(xData')
%     size(yData(k,:))
%     size(xin)
    F = TriScatteredInterp(xData',yData(k,:)');
    out = F(xin');
%     size(yout(k,:))
%     size(out)
    yout(k,:) = out';
end

if (islog)
    yout = logical(round(yout));
end