function qvec = gatherFromSubField(eset, subfield, fieldname, varargin)
%function qvec = gatherField(eset, subfield, fieldname, varargin)
%
%what you would get from qvec = [eset.expt.gatherSubField(subfield, fieldname,
%varargin)] 
%if that were a valid matlab function call
%
%passing 'inds', inds returns only qvec(:,inds) (and inds does not
%propogate to later function calls)
% see Experiment.gatherFromSubField, 
% see Track.getSubFieldDQ and 
%    TrackPart.getDerivedQuantity for optional params

inds = [];
varargin = assignApplicable(varargin);
qvec = [];
if (strcmpi(fieldname, 'expNum'))
    for j = 1:length(eset.expt)
        qv{j} = j*ones(size(eset.expt(j).gatherFromSubField(subfield,'trackNum', varargin{:})));       
    end    
    qvec = [qv{:}];
else
    for j = 1:length(eset.expt)
        qvec = [qvec eset.expt(j).gatherFromSubField(subfield, fieldname, varargin{:})];
    end
end
if (~isempty(inds))
    qvec = qvec(:,inds);
end
