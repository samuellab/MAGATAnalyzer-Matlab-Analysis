function tp = fromFile (tp, fid, loadImage, loadContour, camcalinfo)
%ImTrackPoint.fromFile
%function tp = fromFile (tp, fid, loadImage, loadContour, camcalinfo)
%
%loadContour: ignored
%ts = tic;
intType = 'int32';
%floatType = 'float32';

if (~exist('camcalinfo', 'var'))
    camcalinfo = [];
end
%tic
tp = fromFile@TrackPoint(tp, fid, loadImage, loadContour, camcalinfo);

w = fread(fid,1,intType);
h = fread(fid,1,intType);
if (w ~=0 && h ~= 0)
     tp.imOffset = int16(fread(fid,[1 2],intType));
     if (loadImage)
         tp.imData = transpose(uint8(fread(fid,[w h], 'uchar')));
     else
         fread(fid, [w h], 'uchar');
%         fseek(fid,w*h,'cof'); %move ahead the size of the image
     end
end
%toc(ts)
%pause;
%toc
