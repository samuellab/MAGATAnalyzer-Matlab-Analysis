function [xc, np, tx, nt] = crosscorrelate (expt, fieldname1, fieldname2, varargin)
%cross correlates track fields
%function [xc, np, tx, nt] = crosscorrelate (expt, fieldname1, fieldname2, varargin)
%
%cross correlates track fields (see Track.crosscorrelate)
%
%assumes same interpolation time (dr.interptime) for all tracks
%
%outputs:
%XC is the unnormalized cross correlation (normalized is xc./np)
%TX is the time axis for the cross correlation
%NP is the number of points contributing to a certain bin
%NT is the number of tracks contributing to a certain bin
%
%xc(T) = <dot(a(t),b(t-T))>; -N<T<N and <> denotes the average
%XC is the return 1x(2N-1) vector XC(j) = xc(j-N);
%
%inputs:
%EXPT: a member of the Experiment class
%FIELDNAME1, FIELDNAME2: the names of the fields to cross-correlate
%VARARGIN: any parameter/value pair passed to Track/crosscorrelate, see
%below
%
%we define the cross-correlation to be
%XC(T) = <dot(a(t),b(t-T))>; -N<T<N and <> denotes the average
%xc is the return 1x(2N-1) vector xc(j) = XC(j-N);
%
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


expt.calculateDerivedQuantity('eti');
npts = expt.evaluateTrackExpression('length(track.dq.eti)');
N = max(npts);

xc = zeros([1 2*N-1]);
np = xc;
nt = np;
tx = ((1:(2*N-1)) - N)*expt.track(1).dr.interpTime;

for j = 1:length(expt.track)
    [x, npt] = expt.track(j).crosscorrelate(fieldname1, fieldname2, varargin{:});
    n = (length(x) + 1) / 2;
    if (isempty(x))
        continue;
    end
    xc((N-n+1):(N+n-1)) = xc((N-n+1):(N+n-1)) + x;
    np((N-n+1):(N+n-1)) = np((N-n+1):(N+n-1)) + npt;
    nt((N-n+1):(N+n-1)) = nt((N-n+1):(N+n-1)) + 1;
end

ind = find(np ~= 0);
xc = xc(ind);
np = np(ind);
tx = tx(ind);
nt = nt(ind);