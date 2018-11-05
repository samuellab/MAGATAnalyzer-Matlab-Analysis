function loadTracksFromDir(expt, basedir, varargin)
% function loadTracksFromDir(expt, basedir, varargin)
%
% loads tracks that were saved separately, only if expt.track is empty or 
% optional arg 'overwrite',true or 'append',true is passed
%
% basedir is used only if expt.savedTrackDir does not exist
overwrite = false;
append = false;
existsAndDefault('basedir', '');
varargin = assignApplicable(varargin);

if (~isempty(expt.track) && ~(append || overwrite))
    return;
end

if isdir(fixFileNameWin(expt.savedTrackDir))
    trackdir = expt.savedTrackDir;
else
    if ~isempty(basedir)
        trackdir = fullfile(basedir, expt.savedTrackRelDir);
    else
        return;
    end
end

dd = dir(fullfile(fixFileNameWin(trackdir), '*.mat'));

if (isempty(dd))
    return;
end

for j = 1:length(dd)
    temp = Track.fromMatFile(fixFileNameWin(fullfile(trackdir, dd(j).name)));
    if (isempty(temp))
        continue;
    end
    if (~exist('track', 'var'))
        track = temp;
    else
        try
            track = [track temp]; %#ok<AGROW>
        catch me
            disp (me.getReport());
            error ('track concatenation failed - incompatible types?');
        end
    end
end

if (~exist('track', 'var'))
    return;
end

if (append)
    try %#ok<UNRCH>
        expt.track = [expt.track track];
    catch me
        disp (me.getReport());
        error ('track concatenation failed - incompatible types?');
    end
else
    expt.track = track;
end