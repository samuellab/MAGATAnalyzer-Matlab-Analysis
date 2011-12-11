function tp = fromFile (tp, fid, loadImage, loadContour, camcalinfo)
%MaggotTrackPoint.fromFile
%function tp = fromFile (tp, fid, loadImage, loadContour, camcalinfo)
%

intType = 'int32';
floatType = 'float32';

if (~exist('camcalinfo', 'var'))
    camcalinfo = [];
end
tp = fromFile@ImTrackPoint(tp, fid, loadImage, loadContour, camcalinfo);
tp.targetArea = int16(fread(fid,1,intType));
tp.threshold = single(fread(fid,1,floatType));
tp.htValid = logical(fread(fid,1,intType));
tp.head = single(readPointsFromFile (fid, 1, intType, camcalinfo));
tp.mid = single(readPointsFromFile (fid, 1, intType, camcalinfo));
tp.tail = single(readPointsFromFile (fid, 1, intType, camcalinfo));
ncpts = fread(fid,1,intType);
if (tp.threshold < 0 && tp.targetArea < 0 && ncpts <= 1)
    %ncpts
    %ftell(fid)
    %pause
end
if (loadContour)
    tp.contour = single(readPointsFromFile (fid, ncpts, intType, camcalinfo));
else
    readPointsFromFile (fid, ncpts, intType, camcalinfo);
end

nmpts = fread(fid, 1, intType);
tp.spine = single(readPointsFromFile (fid, nmpts, floatType, camcalinfo));