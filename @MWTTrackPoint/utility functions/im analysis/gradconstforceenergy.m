function gconst = gradconstforceenergy(ctr, closed)
%function ubend = bendenergy(ctr, closed)
%closed if we should connect contour(:,1) with contour(:,end)
existsAndDefault('closed', true);
if (closed)
    dxr = diff(ctr(:, [end 1:end]), [], 2); % (x_j - x_(j-1))
    dxf = diff(ctr(:, [1:end 1]), [], 2); % x_j+1 - x_j
else
    dxr = diff(ctr(:, [1 1:end]), [], 2); % (x_j - x_(j-1))
    dxf = diff(ctr(:, [1:end end]), [], 2); % x_j+1 - x_j
end
%dx = deriv(ctr, sigma, 'padtype', 'circular');

dxfl = sqrt(sum(dxf.^2));
dxrl = sqrt(sum(dxr.^2));
dxfl(dxfl == 0) = 1;
dxrl(dxrl == 0) = 1;

gconst = -dxf./[dxfl;dxfl]+dxr./[dxrl;dxrl];