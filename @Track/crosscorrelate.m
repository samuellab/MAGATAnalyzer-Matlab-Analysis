function [xc, np, tx] = crosscorrelate (track, fieldname1, fieldname2, varargin)
%cross correlates fieldname1 and fieldname2
%function [xc, np, tx] = crosscorrelate (track, fieldname1, fieldname2, varargin)
%
%returns the cross correlation for tracks.dq.(fieldname1) tracks.dq.(fieldname2) 
%we define the cross-correlation to be
%xc(T) = <dot(a(t),b(t-T))>; -N<T<N and <> denotes the average
%XC is the return 1x(2N-1) vector XC(j) = xc(j-N);
%
%outputs:
%XC is the unnormalized cross correlation (normalized is XC./NP)
%NP is the number of points contributing to a certain bin
%TX is the time axis for the cross correlation, so if TX(j) = tau, then
%   XC(j) = xc(tau)
%
%inputs:
%TRACK: a member of the Track class
%FIELDNAME1, FIELDNAME2: the names of the fields to cross-correlate
%VARARGIN: any parameter/value pair below
%
%'row', row number(s)
%   if row = m, just finds the correlation of track.dq.(fieldnames)(m,:)
%'inRuns', true/false
%   if inRuns is true, we take the autocorrelation over the whole track, but
%   first interpolate the fields over only the run indices
%   this is useful if the fields are ill-defined between runs (e.g. velocity
%   direction)
%'withinRuns', true/false
%   if withinRuns is true, we find the correlation only within each run
%'isangle', true/false
%   if isangle is true, then we compute the correlation as cos (theta1 -
%   theta2) instead of as theta1*theta2
%'timerange', [mintime maxtime] -- only consider data from this time range

q1 = track.getDerivedQuantity(fieldname1);
q2 = track.getDerivedQuantity(fieldname2);
eti = track.getDerivedQuantity('eti');
timerange = [];%[min(eti) max(eti)];
isangle = false;    

row = 1:size(q1, 1);
inRuns = false;
withinRuns = false;
varargin = assignApplicable(varargin);
if (isempty(timerange))
    timerange = [min(eti) max(eti)];
end
validtime = eti >= min(timerange) & eti <= max(timerange);


if (isangle)
    if ((size(q1,1)) ~= 1)
        error ('Asked to correlate a field as an angle, but passed a vector field, not a scalar field');
    end
    q1 = unwrap(q1);
    q2 = unwrap(q2);
end
    

if (inRuns)
    if isempty(track.run)
        inds = [];
    else
        inds = find(track.isrun);
    end
    if (isempty(inds))
        xc = [];
        tx = [];
        np = 0;
%        warning('GERSHOW:XC02', 'you asked for values in runs, but track.run is empty or has no points');
        return;
    end
    allinds = min(inds):max(inds);
    extrapval = zeros(size(q1(:,1)));
    if (isangle)
        q1 = interp1(inds, q1(:,inds)', allinds, 'linear');
        q2 = interp1(inds, q2(:,inds)', allinds, 'linear');
    else 
        q1 = interp1(inds, q1(:,inds)', allinds, 'linear','extrap',extrapval)';
        q2 = interp1(inds, q2(:,inds)', allinds, 'linear','extrap',extrapval)';
    end
    validtime = validtime(1,allinds);
end

if (isangle)
    q1 = [cos(q1);sin(q1)];
    q2 = [cos(q2);sin(q2)];
    row = 1:2;
end

if (withinRuns)
    [xc,np,tx] = crossCorrelateRuns(track, q1, q2, row, validtime);
    return
end

if (any(~isfinite(q1(:))) || any (~isfinite(q2(:))))
    xc = [];
    tx = [];
    np = 0;
    warning('GERSHOW:XC01', 'non-finite value encountered: not computing correlation');
    return;
end

[~, xc, np] = xcorrVec(q1(row,validtime), q2(row,validtime));
tx = ((1:length(xc)) - nnz(validtime))*track.dr.interpTime;

function [xc, np, tx, nr] = crossCorrelateRuns(track, q1, q2, row, validtime)

npr = zeros(size(track.run));

    for j = 1:length(track.run)
        npr(j)  = length(track.run(j).inds);
    end
    N = max(npr);

    xc = zeros([1 2*N-1]);
    np = xc;
    nr = np;
    tx = ((1:(2*N-1)) - N)*track.dr.interpTime;

    for j = 1:length(track.run)
        if any(~validtime(track.run(j).inds))
            continue;
        end
        if (any(any(~isfinite(q1(row,track.run(j).inds)))) || any(any(~isfinite(q2(row,track.run(j).inds)))))
            warning('GERSHOW:XC01', 'non-finite value encountered: not computing correlation');
            continue;
        end
        [~, x, npt] = xcorrVec(q1(row,track.run(j).inds), q2(row,track.run(j).inds));
        n = (length(x) + 1) / 2;
        if (isempty(x))
            continue;
        end
    
        xc((N-n+1):(N+n-1)) = xc((N-n+1):(N+n-1)) + x;
        np((N-n+1):(N+n-1)) = np((N-n+1):(N+n-1)) + npt;
        nr((N-n+1):(N+n-1)) = nr((N-n+1):(N+n-1)) + 1;
    end
