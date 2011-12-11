function inds = indsAtTime(track, elapsedTime)
% finds the indices nearest to the elapsedTime(s) specified;
% function inds = indsAtTime(track, elapsedTime)
%
% finds the indices nearest to the elapsedTime(s) specified;
% returns no value for out of range indices (so size(inds) is not
% necessarily the same as size(elapsedTime) )
%
% outputs:
%   INDS: a list of inds s.t. track.dq.eti(inds) approx ELAPSEDTIME
% inputs:
%   TRACK: a member of the Track class
%   ELAPSEDTIME: the times to find indices at

x = track.getDerivedQuantity('eti');
y = 1:length(x);
inds = interp1(x, y, elapsedTime, 'nearest', NaN);
inds = inds(isfinite(inds));