function gspr = gradspringenergy(ctr, closed)
%function ubend = bendenergy(ctr, closed)
%closed if we should connect contour(:,1) with contour(:,end)
existsAndDefault('closed', true);
if (closed)
    dxr = diff(ctr(:, [end 1:end]), [], 2); % (x_j - x_(j-1))
    dxf = diff(ctr(:, [1:end 1]), [], 2); % x_j+1 - x_j
    dx = 0.5 * (dxr + dxf);
   % dx = 0.5 * (diff(ctr(:, [end 1:end]),[],2) + diff(ctr(:, [1:end 1]),[],2));
    ddx = diff(ctr(:, [end 1:end 1]), 2, 2);
else
    dxr = diff(ctr(:, [1 1:end]), [], 2); % (x_j - x_(j-1))
    dxf = diff(ctr(:, [1:end end]), [], 2); % x_j+1 - x_j
    dx = 0.5 * (dxr + dxf);
 %   dx = 0.5 * (diff(ctr(:, [1 1:end]),[],2) + diff(ctr(:, [1:end end]),[],2));
    dx(:,1) = diff(ctr(:,[1 2]),[],2);
    dx(:,end) = diff(ctr(:, [end-1 end]),[],2);
    ddx = 0.5 * (diff(dx(:, [1 1:end]),[],2) + diff(dx(:, [1:end end]),[],2));

end
%dx = deriv(ctr, sigma, 'padtype', 'circular');
dxfl3 = (sum(dxf.^2)).^(1.5);
dxrl3 = (sum(dxr.^2)).^(1.5);
gspr = -ddx;% * length(ctr);% + dxf./[dxfl3;dxfl3]-dxr./[dxrl3;dxrl3];