function str = toMWTTrackString (track, varargin) %writes the appropriate string for a .blobs file
%function str = toMWTTrackString (track, varargin) 
%writes the appropriate string for a .blobs file
%
%track < MaggotTrack

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
% if (~all(valid) || track.nt > 1)    
%     pt = pt);
% end

str = cell(length(track.pt)+1,1);
str{1} = sprintf ('%% %d\n', track.trackNum);
pt = track.pt;
if(isempty(track.expt))
    cc = [];
else
    cc = track.expt.camcalinfo;
end
for j = 1:length(pt)
    str{j+1} = sprintf('%s\n', pt(j).toMWTBlobLine(cc));
end
str = [str{:}];


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
