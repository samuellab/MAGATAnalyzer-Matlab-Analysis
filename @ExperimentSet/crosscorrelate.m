function [xc, np, tx, nt] = crosscorrelate (eset, fieldname1, fieldname2, varargin)
%cross correlates track fields
%function [xc, np, tx, nt] = crosscorrelate (eset, fieldname1, fieldname2, varargin)
%
%cross correlates track fields (see Track.crosscorrelate)
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


maxTime = zeros([1 length(eset.expt)]);
dt = zeros([1 length(eset.expt)]);
for j = 1:length(eset.expt)
    [xc, np, tx, nt] = eset.expt(j).crosscorrelate (fieldname1, fieldname2, varargin{:});
    crosscorr{j} = xc; %#ok<*AGROW>
    numpts{j} = np;
    timeaxis{j} = tx;
    numtracks{j} = nt;
    maxTime(j) = max(tx(np > 0));
    dt(j) = diff(tx(1:2));
end
tx = 0:min(dt):maxTime;
tx = [-tx(end:-1:2) tx]; %this way, tx always includes 0
%tx = -maxTime:(min(dt)):maxTime;
xc = zeros(size(tx));
np = xc;
nt = xc;

for j = 1:length(eset.expt)
    xc = xc + interp1(timeaxis{j}, crosscorr{j}, tx, 'linear', 0);
    np = np + interp1(timeaxis{j}, numpts{j}, tx, 'linear', 0);
    nt = nt + interp1(timeaxis{j}, numtracks{j}, tx, 'linear', 0);
end