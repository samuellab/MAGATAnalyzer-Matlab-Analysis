function camcalinfo = cameraCalibration(realPoints, imagePoints, centerPoint, camcalinfoStart, z)
%function camcalinfo = cameraCalibration(realPoints, imagePoints, centerPoint, camcalinfoStart, z)
%
%REALPOINTS is either a 2xN array of points representing the location of
%corresponding imagePoints or
%a 3xN array of points representing the x,y,z locations of image points or
%a 4x1 (or 1x4) array [nrows ncols rowspacing colspacing]
%
%in the second case, we arrange the imagePoints into an array of
%nrowsxncols
%then we set the real value of the lower left corner to 0,0 and every other
%point is at x = i*colspacing; y = j*rowspacing
%
%
%IMAGEPOINTS is a 2xN or 2x(nrows*ncols) array
%
%centerPoint is the x,y location of the central pixel of the camera CCD 
%our 5 MP cameras are 2592x1944, so the central pixel is 1296, 972
%
%
%CAMCALINFO contains CM, the camera matrix, and K, a measure of radial
%distortion
%CAMCALINFOSTART is an optional parameter that gives a starting point for
%iteration;  only camcalinfoStart.K is used.
%
%Z optional param: defaults to 0;  z location of coplanar real points

if (size(imagePoints,1) ~= 2)
    imagePoints = imagePoints';
end
existsAndDefault('z', 0);

xi = imagePoints(1,:) - centerPoint(1);
yi = imagePoints(2,:) - centerPoint(2);

if (length(realPoints) == 4)
    nrows = realPoints(1);
    ncols = realPoints(2);
    rs = realPoints(3);
    cs = realPoints(4);
    [yi,I] = sort(yi);
    xi = xi(I);
    for j = 1:nrows
        inds = (j-1)*ncols + (1:ncols);
        [x,I] = sort(xi(inds));
        rowi(j).x = x;
        rowi(j).y = yi(inds(I));
        rowr(j).x = (0:(ncols-1)) * cs;
        rowr(j).y = repmat((j-1)*rs,size(x));
    end
    xi = [rowi.x];
    yi = [rowi.y];
    xr = [rowr.x];
    yr = [rowr.y];
    z = repmat(z,size(xr));
  
    
else
    if (size(realPoints,1) ~= 2 && size(realPoints,1) ~= 3)
        realPoints = realPoints';
    end
    xr = realPoints(1,:);
    yr = realPoints(2,:);
    if (size(realPoints,1) == 3)
        z = realPoints(3,:);
    else
        z = repmat(z, size(xr));
    end
end

cmd.K = [0 0 0];
existsAndDefault('camcalinfoStart', cmd);

[CM,K] = cameraMatrixPlusDistortion(xi, yi, xr, yr, z, camcalinfoStart.K);

camcalinfo.CM = CM;
camcalinfo.K = K;
camcalinfo.xc = centerPoint(1);
camcalinfo.yc = centerPoint(2);
