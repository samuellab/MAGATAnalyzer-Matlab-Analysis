function pruneTracks(expt, starttimerange, startrect)
%delete tracks that start outside designated ranges
%function pruneTracks(expt, starttimerange, startrect)
%
%outputs: none
%inputs:
%EXPT: a member of the Experiment class
%STARTTIMERANGE: [t0 t1], the range of valid start times
%STARTRECT: [x0 y0 x1 y1], the space of valid start positions
%pruneTracks removes completely any track that starts outside
%STARTTIMERANGE in elapsedTime, as well as any track that starts
%outside STARTRECT; leave STARTTIMERANGE or STARTRECT empty to disable
    
starttime = transpose(expt.elapsedTime([expt.track.startFrame]+1));

startloc = zeros([2 length(expt.track)]);
for j = 1:length(expt.track)
    startloc(:,j) = expt.track(j).pt(1).loc;
end

valid = true(size(expt.track));
if (exist('starttimerange','var') && ~isempty(starttimerange))
    valid = starttime >= min(starttimerange) & starttime <= max(starttimerange);
end

if (exist('startrect', 'var') && ~isempty(startrect))
    valid = valid & startloc(1,:) >= startrect(1) & startloc(1,:) <= startrect(3) & ...
                startloc(2,:) >= startrect(2) & startloc(2,:) <= startrect(4);
end

expt.track = expt.track(valid);

    