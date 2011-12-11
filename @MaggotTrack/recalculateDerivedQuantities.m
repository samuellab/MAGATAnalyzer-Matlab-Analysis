function recalculateDerivedQuantities(track, varargin)
% clears and recalculates all derived quantities
% function recalculateDerivedQuantities(track, varargin)
% outputs: none
% inputs:
%   TRACK < MaggotTrack
%   VARARGIN: if a list of field names is passed, only those fields are cleared
%       and recalculated

recalculateDerivedQuantities@Track(track,varargin{:});