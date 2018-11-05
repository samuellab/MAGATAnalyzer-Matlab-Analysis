function shapeModel = getTypicalMaggotShape(track)%, varargin)
%function [spineim, rel_width_vs_length, shapex, shapey] = getTypicalMaggotShape(track)

maxpts = 500;
usecamcalinfo = false;
%varargin = assignApplicable(varargin);


if (isempty(track.pt(1).imData))
    disp ('reloading track to get images');
    track.expt.reloadTrack(track.trackNum, 'usecamcalinfo', usecamcalinfo);
    disp ('done');
end

%find maggot track points where body bend angle is small (< 7 deg)
iinds = find(logical(track.getDerivedQuantity('ihtValid')) & abs(track.getDerivedQuantity('spineTheta')) < deg2rad(7));
inds = track.getDerivedQuantity('mapinterpedtopts', false,iinds);
inds = unique(inds);

%if there are too many such points, select randomly from them
if (length(inds) > maxpts)
    inds = inds(randperm(length(inds)) < maxpts);
end


thresh = mean([track.pt(inds).threshold]);
nsp = 100; %number of points along spine for sampling
nwp = 51; %number of points across body for sampling -- changed from 31 to 51 by mhg on 6/19

%mean length of spine (in cm, if calibrated) across sample images
meanSpineLength = mean(track.getDerivedQuantity('spineLength',false,iinds));

%area in pixels / length (in cm, if calibrated)
width_in_pixels = 5*mean(track.getDerivedQuantity('iarea',false,iinds))./meanSpineLength; %changed from 2 to 5 by mhg on 6/25
if (~isempty(track.expt) && usecamcalinfo)
    cc = track.expt.camcalinfo;
else
    cc = [];
end

%if calibrated, correct width pixels^2/cm * cm/pixels
%correct length cm * pixels/cm
if(~isempty(cc))
    width_in_pixels = width_in_pixels/cc.pixelsPerRealUnit;
    meanSpineLength = meanSpineLength*cc.pixelsPerRealUnit;
end



sl = mean(track.getDerivedQuantity('spineLength',false,iinds));
shapey = linspace(-.25*sl, .25*sl,nsp); %change from -0.1,1.1 to -0.25,1.25 by mhg on 6/25, also updated imSpineSampled
pixLen = shapey(end);

%possible bug -- shapey doesn't make a whole lot of sense here
%Spine length ought to be in the proper units already
%if(~isempty(cc))
%    shapey = shapey/cc.pixelsPerRealUnit;
%end

shapex = linspace(-width_in_pixels/2, width_in_pixels/2, nwp);

%sample first track point image along spine
spineim = imSpineSampled (track.pt(inds(1)), nsp, nwp, width_in_pixels,cc);

%timer stuff
ts1 = tic;
etold = 0;

for j = inds(2:end)
    %accumulate spine image by adding on sampled points
    spineim = spineim + imSpineSampled (track.pt(j), nsp, nwp, width_in_pixels,cc);
    
    %timer stuff
    et = toc(ts1);
    if (et - etold > 60)
        etold = et;
        disp ([num2str(et, 2) ' s elapsed. ' num2str(find(inds == j,1)) ' / ' num2str(length(inds))]);
    end
end

%normalize by number of images
spineim = spineim/length(inds);

%eliminate left right assymmetry
ssifull = 0.5*(spineim + spineim(:,end:-1:1)); %symmetrized image
ssifull = blurim(ssifull, 1);
ssi = ssifull(:,ceil((size(ssifull,2)+1)/2):end);

%find first and last y-index (along spine) where there is at least one
%pixel above threshold
ind1 = find(any(ssi >= thresh, 2), 1, 'first');
ind2 = find(any(ssi >= thresh, 2), 1, 'last');

inds = ind1:ind2;
L = diff(shapey([ind1 ind2]));
rel_width_vs_length = zeros(size(inds));

%find width vs. length (outline)
for j = 1:length(inds)
    ind = find(ssi(inds(j),:) >= thresh, 1, 'last');
    if (isempty(ind))
        continue;
    end
    i1 = find(ssi(inds(j),:) < ssi(inds(j),ind), 1, 'first');
    if (isempty(i1) || i1 > inds(j))
        i1 = 1;
    end
    i2 = find(ssi(inds(j),(ind+1):end) > ssi(inds(j),ind), 1, 'first');
    if (isempty(i2))
        i2 = size(ssi,2);
    else
        i2 = ind + i2;
    end
    if (i1 == i2)
        i1 = max(i1 -1, 1);
        i2 = i1+1;
    end
    xdat = ssi(inds(j),i1:i2);
    ydat = i1:i2;
    [xdat,I] = unique(xdat); ydat = ydat(I);
    rel_width_vs_length(j) = interp1(xdat, ydat, thresh, 'spline', 'extrap');
end
%shapey = shapey-shapey(max(ind1-1, 1));
%rel_width_vs_length = double([0 rel_width_vs_length*width_in_pixels/nwp/L 0]);
wvl = zeros(size(shapey));
wvl(inds) = double(rel_width_vs_length*width_in_pixels/nwp);
inds = max((inds(1)-1),1):min(inds(end)+1, length(wvl));
shapey = double(shapey-shapey(inds(1)));

%check math -- some of this is wrong!
shapeModel.image = double(spineim); %not symmetrized!
shapeModel.symmetrizedImage = double(ssifull); %double([fliplr(ssi) ssi]);
shapeModel.spine_coord = shapey;
shapeModel.rel_spine_coord = shapey/shapey(inds(end));
shapeModel.width_coord = double(shapex);

shapeModel.length = shapey(inds)-shapey(inds(1));
shapeModel.pixelLen = pixLen;
shapeModel.widthVsLength = wvl(inds);
    
shapeModel.pixelx = double(shapeModel.width_coord);
shapeModel.pixely = double(shapeModel.rel_spine_coord * meanSpineLength);

function spineim = imSpineSampled (mtp, nsp, nwp, widthInPixels, camcalinfo)
%function imSpineSampled (mtp, nsp, nwp, widthInPixels, camcalinfo)

debug = false;

%spline fit to spine
pp = cscvn(double([mtp.tail mtp.spine mtp.head]));
sp = fnval(pp, linspace(pp.breaks(1), pp.breaks(end), nsp));
l = [0 cumsum(sqrt(sum(diff(sp,[],2).^2)))];
sp = interp1(l,sp',linspace(0,l(end),nsp))';

%extend by 25%
l = [0 cumsum(sqrt(sum(diff(sp,[],2).^2)))];
sp = interp1(l,sp',linspace(-0.25*l(end),1.25*l(end),nsp),'linear', 'extrap')'; %change from 5% by MHG 6/25

%re-regularize
l = [0 cumsum(sqrt(sum(diff(sp,[],2).^2)))];
sp = interp1(l,sp',linspace(0,l(end),nsp))';

if (~isempty(camcalinfo))
    sp = camcalinfo.camPtsFromRealPts(sp);
end

%translate spine to align with image; +1 is because MATLAB is stupid(1-indexed)
sp(1,:) = sp(1,:) - double(mtp.imOffset(1)) + 1;
sp(2,:) = sp(2,:) - double(mtp.imOffset(2)) + 1;

%t is tangent, n is normal
t = diff(sp,[],2);
t = 0.5*(t(:,[1 1:end]) + t(:,[1:end end]));
s = sqrt(sum(t.^2));
n = [t(2,:)./s;-t(1,:)./s];

spineim = zeros(nsp, nwp);

%for each point along spine, interpolate along line locally perpendicular
%to spine at that point
for j = 1:nsp
    pts = repmat(sp(:,j), [1 nwp]) + (n(:,j)*linspace(-widthInPixels/2, widthInPixels/2, nwp));
    
    spineim(j,:) = interp2(double(mtp.imData), pts(1,:),pts(2,:),'*linear', 0);
    if (debug)
        subplot(2,1,1); imagesc(double(mtp.imData)); hold on; plot (sp(1,:), sp(2,:), 'y-', pts(1,:), pts(2,:), 'r-'); colormap gray; hold off
        subplot(2,1,2); imagesc(spineim); pause(0.1);
    end
end