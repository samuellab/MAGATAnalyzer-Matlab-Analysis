function [ac, np, tx] = autocorrelate (track, fieldname, varargin)
%returns the auto correlation for track.dq.(fieldnames)
%function [ac, np, tx] = crosscorrelate (track, fieldname, varargin)
%
%
%outputs:
%AC is the unnormalized auto-correlation (normalized is xc./np)
%NP is the number of points contributing to a certain bin
%TX is the time axis for the cross correlation
%
%we define the auto-correlation to be
%ac(T) = <dot(a(t),a(t-T))>; 0<T<N and <> denotes the average
%AC is the return 1xN vector AC(j) = AC(j-1);
%
%inputs:
%TRACK: a member of the Track class
%FIELDNAME: the name of the field to auto-correlate
%VARARGIN: any parameter/value pair passed to Track/crosscorrelate, see
%below
%
%arguments to pass in:
%'row', row number(s)
%if row = m, just finds the correlation of track.dq.(fieldnames)(m,:)\
%'inRuns', true/false
%if inRuns is true, we take the autocorrelation over the whole track, but
%first interpolate the fields over only the run indices
%this is useful if the fields are ill-defined between runs (e.g. velocity
%direction)
%'withinRuns', true/false
%if withinRuns is true, we find the correlation only within each run
%'timerange', [mintime maxtime] -- only consider data from this time range


[xc, np, tx] = track.crosscorrelate(fieldname, fieldname, varargin{:});
if (~isempty(xc))
    %n = ceil(length(xc)/2);
    n = find(tx == 0, 1);
end
ac = xc(:,n:end);
np = np(n:end);
tx = tx(n:end);
