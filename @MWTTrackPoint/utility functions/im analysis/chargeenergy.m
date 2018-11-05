function ucharge = chargeenergy(ctr, closed)
%function ubend = bendenergy(ctr, closed)
%closed if we should connect contour(:,1) with contour(:,end)
%gives a 1/r potential between nearest neighbors
existsAndDefault('closed', true);
if (closed)
    dx = diff(ctr(:, [end 1:end]),[],2);
 %   dxr = diff(ctr(:, [end 1:end]),[],2);
%    dxf = diff(ctr(:, [1:end 1]),[],2);
   
else
    dx = diff(ctr,[],2);
%    dxr = diff(ctr(:, [1 1:end]),[],2);
%    dxf = diff(ctr(:, [1:end end]),[],2);
end
dl = sqrt(sum(dx.^2));
n = length(dl);
%dx = deriv(ctr, sigma, 'padtype', 'circular');
ucharge = sum(1./dl)/n.^2;