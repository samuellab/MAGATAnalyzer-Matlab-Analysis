function track = fromFile(fid, ptType, loadImageByIndex, loadContour, camcalinfo, minpts)
% loads a track from disk
% function track = fromFile(fid, ptType, loadImageByIndex, loadContour,
% camcalinfo)
%
% outputs:
%   TRACK: a track or a maggot track
% inputs:
%   FID: open file ID
%   PTTYPE: a TrackPoint or subclass thereof;  tells how to load points
%   LOADIMAGEBYINDEX: image(j) is loaded if loadImageByIndex(j) is true
%   LOADCONTOUR: if the maggot contour should be stored on loading (only valid
%       for MaggotPoint)
%   CAMCALINFO: optional information to turn pixel coordinates into real
%       coordinates
%   MINPTS: tracks is not loaded if it has fewer than MINPTS, but the fid
%       pointer is advanced to the next track



if (~exist ('loadImageByIndex', 'var'))
    loadImageByIndex = [];
end
if (~exist ('loadContour', 'var') || isempty (loadContour))
    loadContour = true;
end
if (~exist ('camcalinfo', 'var'))   
    camcalinfo = [];
end
existsAndDefault('minpts', 0);

intType = 'int32';
floatType = 'float32';
if (isa(ptType, 'MaggotTrackPoint'))
    track = MaggotTrack();
    if (~isempty(camcalinfo)) %real points instead of camera points
        track.so.stop_speed_cut = 0.01;
        track.so.start_speed_cut = 0.015;
        track.so.curv_cut = 50;
    end
else
    track = WormTrack();
end
track.locInFile = ftell(fid);
nbytes = fread(fid,1,intType);
startOfTrackPosition = ftell(fid);
track.npts = fread(fid, 1, intType);
if (track.npts >= minpts)
    
    if (length (loadImageByIndex) ~= track.npts)
        if (length(loadImageByIndex) == 1)
            loadImageByIndex = repmat(loadImageByIndex, [1 track.npts]);
        else
            if (~isempty(loadImageByIndex))
                disp(['length of loadImageByIndex = ' num2str(length(loadImageByIndex)) ' but track.npts = ' num2str(track.npts)]);
            end
            loadImageByIndex = repmat(false, [1 track.npts]);
        end
    end
    bob = repmat(ptType, [1 track.npts]);
    %ts = tic;

    for j = 1:track.npts
   
    % track.pt(j) = track.pt(j).fromFile(fid, loadImageByIndex(j), loadContour, camcalinfo);    
        bob(j) = ptType.fromFile(fid, loadImageByIndex(j), loadContour, camcalinfo);   
    
    end
    track.pt = bob;
    %toc(ts);
    bytesLoaded = ftell(fid)-startOfTrackPosition;
    if (nbytes ~= bytesLoaded)
        track.locInFile
        startOfTrackPosition
        track.npts
        bob
        disp (['Expected to load ' num2str(nbytes) ' but loaded ' num2str(bytesLoaded)]);
        pause
    end
    track.startFrame = track.pt(1).ind;
    track.endFrame = track.pt(end).ind;
else
    fseek(fid, startOfTrackPosition + nbytes, 'bof');
end