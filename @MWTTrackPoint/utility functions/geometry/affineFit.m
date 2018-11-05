function [xc, yc, xr, yr, A] = affineFit (xc, yc, xr, yr)
%function [xc, yc, xr, yr, A] = affineFit (xc, yc, xr, yr)


midx = mean(xc);
midy = mean(yc);

[~,I] = sort((xc - midx).^2 + (yc - midy).^2); %find nearest points to center

ncpts = 9;

xro = xr(I(1));
yro = yr(I(1));

xxc = xc(I(1:ncpts));% - xc(I(1));
yyc = yc(I(1:ncpts));% - yc(I(1));
xxr = round(xr(I(1:ncpts)) - xro);
yyr = round(yr(I(1:ncpts)) - yro);

cp = [xxc;yyc]'; 
rp = [xxr;yyr]'; 
R = cp2tform(cp,rp, 'projective');


for j = 1:5

    [xxr,yyr] = tformfwd(R, xc, yc);
    inds = [];
    delta = 0.01 * j;
    while (length(inds) < min(ncpts, length(xc) / 10))
        delta = delta * 2;
        inds = find(abs(xxr - round(xxr)) < delta & abs(yyr - round(yyr)) < delta);
    end

    cp = [xc(inds);yc(inds)]';
    rp = round([xr(inds);yr(inds)])';

    R = cp2tform(cp,rp, 'projective');
end

[xxr,yyr] = tformfwd(R, xc, yc);

% delta = 0.4;
% inds = find(abs(xxr - round(xxr)) < delta & abs(yyr - round(yyr)) < delta);
% 
% xc = xc(inds);
% yc = yc(inds);
xr = round(xxr) + round(xro);
yr = round(yyr) + round(yro);

valid = true(size(xr));
%kill duplicate points by keeping the closest
for j = 1:length(xr)
    if nnz(xr == xr(j) & yr == yr(j)) > 1
        inds = find(xr == xr(j) & yr == yr(j));
        valid(inds) = false;
        [~,I] = min(abs(xxr(inds) - round(xxr(inds))) + abs(yyr(inds) - round(yyr(inds))));
        valid(inds(I)) = true;
    end
end

xc = xc(valid);
xr = xr(valid);
yc = yc(valid);
yr = yr(valid);

