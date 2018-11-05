function segmentTrack (track, wormSegmentOptions, varargin)
% segments TRACK into runs & reorientations using WORMSEGMENTOPTIONS
% function segmentTrack (track, wormSegmentOptions)
%
% worm track segmentation
% outputs: none
% inputs:
%   TRACK < Track
%   WORMSEGMENTOPTIONS < WormSegmentOptions
%
% for explanation of segment options, doc WormSegmentOptions
%
% 'UseExistingSharpTurns', [false]/true, if true, don't find new
%   sharp turns -- so that when we load them from file, we don't nuke
%   existing sharp turns -- if false, we overwrite existing turns, unless
%   user has already started flagging
%
% 'OverwriteExistingSharpTurns', [false]/true, if true, we
%   overwrite existing turns EVEN IF USER HAS ALREADY FLAGGED THEM,
%   eliminating user codes

existsAndDefault('wormSegmentOptions', track.so);
wso = wormSegmentOptions;
track.so = wso;
track.run = [];
track.reorientation = [];
track.sharpTurn = [];
track.flagReorientations(wso, varargin{:});

if (isempty(track.reorientation))
    runstart = 1;
    runend = length(track.dq.eti);
else

    %kill first reorientation if it does not follow a run
    if (track.dq.eti(track.reorientation(1).startInd) - track.dq.eti(1) < wso.minRunTime)
        rs0 = track.reorientation(1).endInd + 1;
        track.reorientation = track.reorientation(2:end);
    else
        rs0 = 1;
    end
    
    %kill last reorientation if it does not precede a run
    if (isempty(track.reorientation))
        runstart = 1;
        runend = length(track.dq.eti);
    else
        if (track.dq.eti(end) - track.dq.eti(track.reorientation(end).endInd) < wso.minRunTime)
            relast = track.reorientation(end).startInd - 1;
            track.reorientation = track.reorientation(1:(end-1));
        else
            relast = length(track.dq.eti);
        end
        runstart = [rs0 ([track.reorientation.endInd]+1)];
        runend = [([track.reorientation.startInd]-1) relast];
    end
end
%run(1) has no previous reorientation
%run(1): next reorientation is reorientation(1)
for j = 1:length(runstart)
    run(j) = Run(track, runstart(j), runend(j));
    %{
    if (j > 1)
        run(j).previousRun = run(j-1);
        run(j).prevReorientation = track.reorientation(j-1);
        run(j-1).nextRun = run(j);
        run(j-1).nextReorientation = track.reorientation(j-1);
    end
    %}
end

track.run = run;
track.isrun = false(size(track.dq.eti));
track.isrun([track.run.inds]) = true;
track.reorientation.calculateMetrics;
