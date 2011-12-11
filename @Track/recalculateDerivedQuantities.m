function recalculateDerivedQuantities(track, varargin)
% clears and recalculates all derived quantities
% function recalculateDerivedQuantities(track, varargin)
% outputs: none
% inputs:
%   TRACK < Track
%   VARARGIN: if a list of field names is passed, only those fields are cleared
%       and recalculated

if (isempty(track.dq))
    return;
end

fnames = fieldnames(track.dq);
valid = track.validDQName(fnames);
vf = fnames(valid);
if (~isempty(varargin))
    if (iscell(varargin{1}))
        vf = intersect(vf, varargin{1});
    else
        vf = intersect(vf, varargin);
    end
end

if (isempty(vf))
    return;
end
%vf
track.dq = rmfield(track.dq, vf);
%track.dq = [];

track.calculateDerivedQuantity(vf);
