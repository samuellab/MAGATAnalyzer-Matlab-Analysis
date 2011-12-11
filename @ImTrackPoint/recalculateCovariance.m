function pt = recalculateCovariance(pt, targetArea, fid)
%function pt = recalculateCovariance(pt, targetArea, fid)

if (~exist('targetArea','var'))
    warning ('ITP:RCC', 'need a targetArea to find object');
    return;
end

if (~exist('fid','var') && isempty(pt.imData))
     warning ('ITP:RCC', 'no image data & no reload info');
    return;
end

if (isempty(pt.imData))
    fseek(fid, pt.locInFile,'bof');
    pt2 = pt.fromFile(fid, true, true, []);
    im = pt2.imData;
else
    im = pt.imData;
end

[bwim, thresh] = thresholdToTargetArea(im, targetArea,pt.loc' - double(pt.imOffset) + 1);
im = im - thresh;
im(~bwim) = 0;
[c,com] = imCov(im);
pt.cov = [c(1,1);c(1,2);c(2,2)];
pt.loc = (com + double(pt.imOffset) - 1)';
pt.area = nnz(bwim(:));