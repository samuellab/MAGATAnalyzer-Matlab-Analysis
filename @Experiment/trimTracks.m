function trimTracks(expt, timerange, validrect)
%removes bad parts of tracks
%function trimTracks(expt, timerange, validrect)
%
%inputs:
%EXPT: a member of the Experiment class
%TIMERANGE: [t0 t1], the range of times (seconds)
%VALIDRECT: [x0 y0 x1 y1], the space of valid positions
%
%trimTracks cuts out any part of the track outside
%min(timerange),max(timerange) and removes any part of the track
%from the point the track leaves validrect until the end of the
%track; leave timerange or validrect empty to disable
if (isempty(timerange) && isempty(validrect))
    return;
end
rval = expt.executeTrackFunction('trim', timerange, validrect);
if (iscell(rval))
    todelete = cell2mat(rval);
else
    todelete = rval;
end
expt.track = expt.track(~todelete);
minpts = (expt.track(1).dr.derivTime + expt.track(1).dr.smoothTime) / expt.track(1).dr.interpTime + 1;
expt.track = expt.track([expt.track.npts] > minpts);

expt.executeTrackFunction('recalculateDerivedQuantities');
expt.assignGlobalQuantities();