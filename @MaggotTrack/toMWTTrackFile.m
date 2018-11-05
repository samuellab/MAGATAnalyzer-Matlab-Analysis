function toMWTTrackFile (track, filestub, varargin) %writes a .blob file for 1 track
%function toMWTTrackFile (track, filestub, varargin) %writes a .blob file for 1 track
%
%writes track to an MWT blob file
%track < MaggotTrack
%filestub - everything before the _NNNNN.blob in the blob file name (incl.
%           directory)
%

if (length(track) > 1)
    for j = 1:length(track)
        track(j).toMWTTrackFile(filestub, varargin{:});
    end
    return;
end

% if (track.nt > 1)    
%     pt = fillInBlanks(track.pt);
% else
%     pt = track.pt;
% end
pt = track.pt;
valid = [pt.area] >= 0.4 * [pt.targetArea];
if (mean(valid) < 0.8)
    valid = [pt.area] > 20; %minimum 20 pixels in area
    if (mean(valid) < 0.8)
        warning ('MGT:TOMWTTRACK', ['a large number of points in track ' num2str(track.trackNum) ' fall below threshold area']);
        return;
    end
end
pt = fillInBlanks(pt(valid));
fid = fopen ([filestub sprintf('_%05d.blob', track.trackNum)], 'wt');
if(isempty(track.expt))
    cc = [];
else
    cc = track.expt.camcalinfo;
end
for j = 1:length(pt)
    fprintf(fid, '%s\n', pt(j).toMWTBlobLine(cc));
end
fclose(fid);

function pt = fillInBlanks(pt)

while (any (diff([pt.ind]) > 1))
    [~,I] = find(diff([pt.ind]) > 1, 1, 'first');
    npts = double(pt(I+1).ind-pt(I).ind - 1);
    newpt = repmat(pt(I), 1, npts);
    for j = double(1:npts)
        newpt(j).ind = pt(I).ind + j;
        newpt(j).loc = (pt(I).loc *(npts+1-j)+pt(I+1).loc *(j))/(npts+1);
    end
    pt = [pt(1:I) newpt pt((I+1):end)];
end
