function im = overlayImageOnExisting (pt, imx, imy, im, camcalinfo, varargin)
%function im = overlayImageOnExisting (pt, imx, imy, im, camcalinfo, varargin)
%PT < ImTrackPoint
%IMX - image x axis
%IMY - image y axis
%IM - image -- size = [length(imy), length(imx), nchannels]
%
%optional arguments, with defaults
%fid = [];
%mode = 'add' or 'a'; %if 'replace' or 'r', replace the data in im, 
%                      if 'a', add to it
%                      if 'max', take maximum of im, pt.imData
%channel = []; % if a number, only replace data in that channel
%
%clipToContour = false; if true, restrict image to region in and near the
%                      contour; only active if pt has a contour or 
%                      if contour is passed as an extra argument
%contour = (pt.contour) or [] ; a 2xN contour of the animal

fid = [];
clipToContour = true;
existsAndDefault('camcalinfo', []);
mode = 'a';
if (any(strcmp('contour', properties(pt))))
    contour = double(pt.contour);
else
    contour = [];
end
channel = [];

varargin = assignApplicable(varargin);


if (isempty(pt.imData))
    if (~isempty(fid) && fid > 0)
        try
            fseek(fid, pt.locInFile, -1);
            pt2 = pt.fromFile(fid, true, true, camcalinfo);
            pt.imData = pt2.imData;
        catch me
            warning('ITP:DTI' , me.getReport);
            return;
        end
    else
        return;
    end
end

x = double(pt.imOffset(1)-1) + (1:size(pt.imData,2));
y = double(pt.imOffset(2)-1) + (1:size(pt.imData,1));

if (~isempty(camcalinfo))
   
    rpts = camcalinfo.realPtsFromCamPts([x(1) x(end); y(1) y(end)]);
    x = linspace(rpts(1,1), rpts(1,end), length(x));
    y = linspace(rpts(2,1), rpts(2,end), length(y));
end

uu = find(imx >= min(x), 1, 'first'):find(imx <= max(x), 1, 'last');
vv = find(imy >= min(y), 1, 'first'):find(imy <= max(y), 1, 'last');

if (isempty(uu) || isempty(vv)) %no overlap
    return;
end

[uu,vv] = meshgrid(uu,vv);
imd = interp2(x,y,double(pt.imData), imx(uu(:)), imy(vv(:)), '*linear', NaN);

if (clipToContour && ~isempty(contour))
    pointsOutside = 2;
    mask = reshape(inpolygon(imx(uu(:)), imy(vv(:)), contour(1,:), contour(2,:)), size(uu));
    mask = imdilate(mask, ones(2*pointsOutside+1));
    %mask = zeros(size(uu));
    imd(~mask(:)) = NaN;
end

valid = isfinite(imd(:));
uu = uu(valid);
vv = vv(valid);
imd = imd(valid);

uu = uu(:); vv = vv(:); imd = imd(:);
% 
% 
% 
% xx = interp1(x, 1:length(x), imx(uu(:)), '*nearest', 'NaN');
% yy = interp1(y, 1:length(y), imy(vv(:)), '*nearest', 'NaN');
% 
% 
% [xx,yy] = meshgrid(x,y);
% x = interp1(imx, 1:length(imx), xx(:), '*nearest', NaN);
% y = interp1(imy, 1:length(imy), yy(:), '*nearest', NaN);
% valid = isfinite(x) & isfinite(y);
% imd = double(pt.imData(valid));
% x = x(valid);
% y = y(valid);
if (isempty(channel))
    channel = 1:size(im,3);
end
% 
% for j = 1:length(channel)
%     plane{j} = im(:,:,channel(j));
% end
% xx = xx(:);
% yy = yy(:);
% imd = imd(:);

for j = 1:length(channel)
%     inds = sub2ind(size(im), y, x, repmat(channel(j),size(x)));
    inds = sub2ind(size(im), vv, uu, repmat(channel(j),size(vv)));
    switch(mode)        
        case {'a', 'add'}
            im(inds) = im(inds) + imd;
        case {'r', 'replace'}
            im(inds) = imd;
        case {'m', 'max'}
            im(inds) = max(im(inds), imd);
    end
end


