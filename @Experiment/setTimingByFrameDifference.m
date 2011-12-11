function setTimingByFrameDifference(expt, deltaT, overrideExisting) 
%function setTimingByFrameDifference(expt, deltaT, overrideExisting) 
%
%sets the elapsedTime to be 0:1:numFrames * deltaT
%if there is already an existing elapsedTime, then does not overwrite it,
%unless overrideExisting is true

existsAndDefault('overrideExisting', false);
existsAndDefault('deltaT', 1);
if (~isempty(expt.elapsedTime) && ~overrideExisting)
    return;
end
tr = expt.track;
maxFrame = max([tr.endFrame]);

expt.dr.interpTime = deltaT;
[tr.dr] = deal(expt.dr);

t = double(0:1:maxFrame)*deltaT;
expt.elapsedTime = double(t');
indx = (1:length(t)) - 1;
for j = 1:length(expt.track)
    expt.track(j).addTime (indx, expt.elapsedTime);
end

expt.timfname = 'asserted timing from inter-frame interval';
