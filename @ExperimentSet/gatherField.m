function qvec = gatherField(eset, fieldname, varargin)
% gathers fieldname from all tracks in all experiments
% function qvec = gatherField(eset, fieldname, varargin)
%
% what you would get from qvec = [eset.expt.gatherField(fieldname,
% varargin)] 
% if that were a valid matlab function call
% inputs: 
% ESET: a member of the ExperimentSet class
% FIELDNAME: either a property of the Track/MaggotTrack/etc. in ESET.expt or
%   a fieldname that can be passed to Track/getDerivedQuantity
% VARARGIN: if fieldname is passed to getDerivedQuantity, VARARGIN{:}
%   is also passed, so any valid arguments to getDerivedQuantity can be
%   appended here
%   e.g. add 'runs', 'reorientations', 'headswings','runend' to get field just in
%   runs/reos/hs
%   'validname', fieldname
%   'validoperation', op -- @(x) logical(x)
%       gathers from the subset of points that stasify
%       that satisfy op(gatherField(fieldname,varargin{:})) == true
%   'inds' -- a list of indices to select

inds = [];
varargin = assignApplicable(varargin);
qvec = [];
for j = 1:length(eset.expt)
    qvec = [qvec eset.expt(j).gatherField(fieldname, varargin{:})];
end

if (~isempty(inds))
    qvec = qvec(:,inds);
end