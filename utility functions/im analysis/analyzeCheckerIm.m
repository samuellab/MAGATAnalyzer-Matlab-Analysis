function analyzedIms = analyzeCheckerIm(checkim, bordersize, cornersize)
%function analyzedIms = analyzeCheckerIm(checkim, bordersize, cornersize)

%first we need to create a bw checkerboard -- let's autosegment using
%k-means

existsAndDefault('cornersize', bordersize);

[idx,c] = kmeans(double(checkim(:)), 2, 'start', double([min(checkim(:)); max(checkim(:))]));
if (c(2) > c(1))
    high = 2;
else
    high = 1;
end
analyzedIms.bwim = imclose(reshape(idx == high, size(checkim)),ones(5));

analyzedIms.distToLight = bwdist(analyzedIms.bwim);
analyzedIms.distToDark = bwdist(~analyzedIms.bwim);
analyzedIms.distToBorder = -analyzedIms.distToLight + analyzedIms.distToDark;

[xd,yd] = imgradient(analyzedIms.distToBorder, 3);
%{
mag = sqrt(xd.^2 + yd.^2);
xd = xd./mag;
yd = yd./mag;
%}

analyzedIms.dirToBorder = atan2(yd,xd);

[cim, cpts] = myharris(analyzedIms.bwim,3);
analyzedIms.cornerMetric = cim;
analyzedIms.cornerLocations = cpts;

cpts = max(round(cpts), 1);
cpts(1,:) = min(cpts(1,:), size(checkim, 2));
cpts(2,:) = min(cpts(2,:), size(checkim, 1));

cim = false(size(checkim));
tinds = sub2ind(size(checkim), cpts(2,:), cpts(1,:));
cim(tinds) = true;
onborder = abs(analyzedIms.distToBorder) <= 1;


%detect lines to determine most likely border orientations
[h,theta] = hough(onborder);

%by squaring the hough transform, we emphasize points of high concentration
[~,I] = sort(sum(h.^2), 'descend');
theta = theta(I);
thetam = theta(1);

%make structuring elements to dilate corner region in an X - shape

semaj1 = strel('line', 2*cornersize, thetam + 45);
semin1 = strel('line', ceil(cornersize/2), thetam - 45);

semaj2 = strel('line', 2*cornersize, thetam - 45);
semin2 = strel('line', ceil(cornersize/2), thetam + 45);

cim1 = imdilate(imdilate(cim, semaj1), semin1);
cim2 = imdilate(imdilate(cim, semaj2), semin2);

cim = imclose(cim1 | cim2, ones(round(cornersize/2)));

%use morphological operations to find border & corner without assumptions
%about orientation of checkerboard

%se = strel('disk', ceil(bordersize/2), 0);

%onborder = imclose(imdilate(onborder,se),ones(5));
%onborder = imerode(imdilate(onborder,se),

onborder = abs(blurim(analyzedIms.distToBorder,3)) < bordersize/2;

analyzedIms.border = onborder &~ cim;
analyzedIms.oncorner = cim & onborder;
analyzedIms.interior = ~(analyzedIms.border | analyzedIms.oncorner);


end

