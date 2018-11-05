function [midpt, dirvec, squareError, fitpts] = fitLine (pts)
%function [midpt, dirvec, squareError, fitpts] = fitLine (pts)
%pts is a DxNPTSxNMEAS; midpt, dirvec are DxNMEAS, squareError is
%1xNMEAS
%(D = dimension, often 2 or 3);  fits a line for each measurement;  assumes
%that pts are evenly spaced along each line

midpt = mean(pts, 2);
po = pts - repmat(midpt, [1 size(pts,2) 1]);

s = repmat(1:size(pts,2), [size(pts,1), 1, size(pts,3)]);
s = s - repmat(mean(s,2), [1, size(pts,2), 1]);
ss = sum(s(1,:,1).^2);
dir = sum(s.*po, 2)/ss;
pf = s.*repmat(dir, [1, size(pts,2), 1]);
fitpts = pf + repmat(midpt, [1 size(pts,2) 1]);
squareError = sum(squeeze(sum((po-pf).^2)));

dir = squeeze(dir);
dirlen = sqrt(sum(dir.^2));
dirvec = dir./repmat(dirlen, [size(dir,1), 1]);

midpt = squeeze(midpt);

%{

perpvec = zeros(size(po));

dirvec = zeros(size(midpt));

for j = 1:size(pts, 3)
    c = cov(po(:,:,j)');
    [v,d] = eig(c);
    %I can't find anywhere in the documentation that says the eigenvalues
    %are sorted, so don't take a chance
    if (d(1,1) > d(2,2))
        dv = v(:,1);
        pv = v(:,2);
    else
        dv = v(:,2);
        pv = v(:,1);
    end
    dirvec(:,j) = dv;
    perpvec(:,:,j) = repmat(pv, [1 size(perpvec,2)]);
end

dist = squeeze(sum(po.*perpvec));
squareError = sum(dist.^2);
 %}