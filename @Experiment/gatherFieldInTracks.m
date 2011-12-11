function qvec = gatherFieldInTracks(expt, fieldname, tracknums, varargin)
%gets all values of 'fieldname' for a subset of tracks in expt
%function qvec = gatherFieldInTracks(expt, fieldname, tracknums, varargin)
%outputs: 
%QVEC: a kxN array of values, where N is the total number of points
%   and k is the dimension of the values of fieldname
%
%inputs: 
%EXPT: a member of the Experiment class
%FIELDNAME: either a property of the Track/MaggotTrack/etc. in EXPT or
%   a fieldname that can be passed to Track/getDerivedQuantity
%TRACKNUMS: the indices of tracks from which to gather field
%VARARGIN: if fieldname is passed to getDerivedQuantity, VARARGIN{:}
%   is also passed, so any valid arguments to getDerivedQuantity can be
%   appended here
%   e.g. add 'runs', 'reorientations', 'headswings' to get field just in
%   runs/reos/hs

if (any(strcmp(fieldname, fieldnames(expt.track(1)))))
    qvec = [expt.track(tracknums).(fieldname)];
else 
    qvec = [];
    for j = tracknums;
        qvec = [qvec expt.track(j).getDerivedQuantity(fieldname, false, varargin{:})];
        if (isempty(qvec))
            break;
        end
    end
end