function qvec = gatherSubField (eset, field, subfield, varargin)
%function qvec = gatherSubField (eset, field, subfield, varargin)
%
%what you would get from [eset.expt.gatherSubField(field,subfield,varargin{:})] if 
%that were valid syntax
%
%outputs:
%QVEC: a kxN array where k is the dimension of field.subfield and N is the
%   number of points
%
%inputs:
%ESET: a member of the ExperimentSet class
%FIELD: a property of ESET.expt.track
%   if field is 'firsths', gathers the first headsweep from every
%   reorientation only
%SUBFIELD: a property of ESET.expt.track.FIELD
%   i.e. if FIELD is 'run', subfield might be 'startTheta', in which case
%   the starting angle of every run is returned
%optional inputs:
%'expandToInds', true/false
%if you pass 'expandToInds',true the subfield is expanded so that
%the return vector is the same length as gatherFromSubField would return
%field must have a subfield inds
%
%for example say a field has a subfield "fieldnumber" = 1 for first field, 2
%for second and so on
%and field(1).inds = [1 2 3], field(2).inds = [1 2 3 4 5], field(3).inds =
%[187]
%
%then gatherSubField(field, fieldnumber, 'expandToInds', true) would yield
%[1 1 1 2 2 2 2 2 3]

qvec = [];
for j = 1:length(eset.expt)
    qvec = [qvec eset.expt(j).gatherSubField(field, subfield, varargin{:})];
end
