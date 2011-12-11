function uspr = springenergy(ctr, closed)
%function ubend = bendenergy(ctr, closed)
%closed if we should connect contour(:,1) with contour(:,end)
existsAndDefault('closed', true);
if (closed)
    dx = 0.5 * (diff(ctr(:, [end 1:end]),[],2) + diff(ctr(:, [1:end 1]),[],2));
   
else
    dx = 0.5 * (diff(ctr(:, [1 1:end]),[],2) + diff(ctr(:, [1:end end]),[],2));
    dx(:,1) = diff(ctr(:,[1 2]),[],2);
    dx(:,end) = diff(ctr(:, [end-1 end]),[],2);
 end
%dx = deriv(ctr, sigma, 'padtype', 'circular');
dL = mean(sqrt(sum(dx.^2)));
A = dL.^2; B = 1/dL;
%uspr = 0.5 * sum(sum(dx.^2)) * length(ctr);% + sum(1./sqrt((sum(dx.^2))));
%uspr = 0.5*A*sum(sum(dx.^2)) + B*sum(1./sqrt((sum(dx.^2)))) - length(ctr);
uspr = 0.5*sum(sum(dx.^2));