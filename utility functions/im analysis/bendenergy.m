function ubend = bendenergy(ctr, closed)
%function ubend = bendenergy(ctr, closed)
%closed if we should connect contour(:,1) with contour(:,end)
existsAndDefault('closed', true);
if (closed)
    dx = 0.5 * (diff(ctr(:, [end 1:end]),[],2) + diff(ctr(:, [1:end 1]),[],2));
    ddx = diff(ctr(:, [end 1:end 1]), 2, 2);
    dl = sqrt(sum(dx.^2));
    that = dx./[dl;dl];
    dthat = 0.5 * (diff(that(:, [end 1:end]),[],2) + diff(that(:, [1:end 1]),[],2));
else
    dx = 0.5 * (diff(ctr(:, [1 1:end]),[],2) + diff(ctr(:, [1:end end]),[],2));
    dx(:,1) = diff(ctr(:,[1 2]),[],2);
    dx(:,end) = diff(ctr(:, [end-1 end]),[],2);
    ddx = [[0;0] diff(ctr, 2, 2) [0;0]];
    dl = sqrt(sum(dx.^2));
    that = dx./[dl;dl];
    dthat = 0.5 * (diff(that(:, [1 1:end]),[],2) + diff(that(:, [1:end end]),[],2));
    dthat(:,1) = diff(that(:,[1 2]),[],2);
    dthat(:,end) = diff(that(:, [end-1 end]),[],2);
    
end
%dx = deriv(ctr, sigma, 'padtype', 'circular');
dl = sqrt(sum(dx.^2));
ds = dl/sum(dl);
%figure(1); 
%plot (1:length(ctr), sum(ddx.^2)/max(sum(ddx.^2)), 1:length(ctr), dl./max(dl))
%cs = sum((ddx.^2))./dl.^3;%./(dl.^2));
%figure(2);
%plotColorLine (ctr(1,:), ctr(2,:), cs, jet);
cs = sum((dthat./[dl;dl])).^2;

ubend = 0.5*sum(cs.*ds);
%ubend = 0.5*sum(sum(ddx.^2));