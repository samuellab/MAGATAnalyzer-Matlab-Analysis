function [ac, np, tx, nt] = autocorrelate (expt, fieldname, varargin)
%finds the autocorrelation of a field with itself
%function [ac, np, tx, nt] = autocorrelate (expt, fieldname, varargin)
%
%
%auto correlates track field (see Track.crosscorrelate)
%
%outputs:
%AC is the unnormalized auto-correlation (normalized is xc./np)
%TX is the time axis for the cross correlation
%NP is the number of points contributing to a certain bin
%NT is the number of tracks contributing to a certain bin
%
%we define the auto-correlation to be
%ac(T) = <dot(a(t),a(t-T))>; 0<T<N and <> denotes the average
%AC is the return 1xN vector AC(j) = AC(j-1);
%
%inputs:
%EXPT: a member of the Experiment class
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
%'isangle', true/false
%   if isangle is true, then we compute the correlation as cos (theta1 -
%   theta2) instead of as theta1*theta2
%'timerange', [mintime maxtime] -- only consider data from this time range


[xc,np,tx,nt] = expt.crosscorrelate(fieldname, fieldname, varargin{:});
N = (length(xc) + 1) / 2;
ac = xc(N:end);
np = np(N:end);
tx = tx(N:end);
nt = nt(N:end);