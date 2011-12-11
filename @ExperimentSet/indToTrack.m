function varargout = indToTrack(eset, ind)
% converts a position in list of all tracks into a track or expt, track indices
% function varargout = indToTrack(eset, ind)
%
% [exptind,tracknum] = eset.indToTrack(ind)
% trackptr = eset.indToTrack(ind)
%
% say t = [eset.expt.track];
% then trackptr = t(ind) = eset.expt(exptind).track(ind)
%
% outputs:
% if 1 output: 
%    TRACKPTR = indicated track (handle)
% if 2 outputs:
%   EXPTIND, TRACKNUM: track is found at ESET.expt(EXPTIND).TRACKNUM
%
% inputs:
% ESET: a member of the ExperimentSet class
% IND: the index of the track in [eset.expt.track]

if (nargout == 1)
    t = [eset.expt.track];
    varargout{1} = t(ind);
    return
end

eind = zeros(size([eset.expt.track]));
tind = eind;
k = 0;
for j = 1:length(eset.expt)
    tinds = 1:length(eset.expt(j).track);
    eind(k + tinds) = j;
    tind(k + tinds) = tinds;
    k = k + length(eset.expt(j).track);
end

varargout{1} = eind(ind);
varargout{2} = tind(ind);
