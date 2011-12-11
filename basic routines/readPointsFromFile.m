function pts = readPointsFromFile (fid, npts, type, camcalinfo)
%function pts = readPointsFromFile (fid, npts, type, camcalinfo)
%
%reads in a set of points (2 x npts of type) from a file then transforms
%them into real coordinates using camcalinfo
%npts
%type
if (npts < 1)
    pts = [];
    return
end
if (npts > 1)
    pts = fread (fid,[2 double(npts)], type);
else
    pts = fread(fid, 2, type);
end
if ~isempty(camcalinfo)
    if (isa (camcalinfo, 'CameraCalibration'))
        pts = camcalinfo.realPtsFromCamPts(pts);
    else
        for j = 1:npts
            [xx,yy] = realPointFromCameraPoint(pts(1,j)-camcalinfo.xc,pts(2,j)-camcalinfo.yc,camcalinfo.CM, camcalinfo.K);
            pts(1,j) = xx;
            pts(2,j) = yy;
        end
    end
end
%{
for j = 1:npts
    if (~isempty(camcalinfo)) 
        [xx,yy] = realPointFromCameraPoint(pts(1,j)-camcalinfo.xc,pts(2,j)-camcalinfo.yc,camcalinfo.CM, camcalinfo.K);
    else
        xx = pts(1,j);
        yy = pts(2,j);
    end
    pts(1,j) = xx;
    pts(2,j) = yy;
end
%}