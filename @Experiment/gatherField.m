function qvec = gatherField(expt, fieldname, varargin)
%gets all values of 'fieldname' for all tracks in expt
%function qvec = gatherField(expt, fieldname, varargin)
%outputs: 
%QVEC: a kxN array of values, where N is the total number of points
%   and k is the dimension of the values of fieldname
%
%inputs: 
%EXPT: a member of the Experiment class
%FIELDNAME: either a property of the Track/MaggotTrack/etc. in EXPT or
%   a fieldname that can be passed to Track/getDerivedQuantity
%VARARGIN: if fieldname is passed to getDerivedQuantity, VARARGIN{:}
%   is also passed, so any valid arguments to getDerivedQuantity can be
%   appended here
%   e.g. add 'runs', 'reorientations', 'headswings' to get field just in
%   runs/reos/hs
%   'validname', fieldname
%   'validoperation', op -- @(x) logical(setNonFiniteToZero(x))
%       gathers from the subset of points that stasify
%       that satisfy op(gatherField(fieldname,varargin{:})) == true

validname = [];
expandToInterped = false;
validoperation = @(x) logical(setNonFiniteToZero(x));
varargin = assignApplicable(varargin);
if (any(strcmp(fieldname, fieldnames(expt.track(1)))))
    if (expandToInterped && all(size(expt.track(1).(fieldname)) == 1))
        qvec = [];       
        for j = 1:length(expt.track)
            qvec = [qvec expt.track(j).(fieldname)*ones(size(expt.track(j).getDerivedQuantity('eti', false, varargin{:})))];
        end
        
    else
        qvec = [expt.track.(fieldname)];
    end
else 
    qvec = [];    
    
    for j = 1:length(expt.track)
        
        qvec = [qvec expt.track(j).getDerivedQuantity(fieldname, false, varargin{:})];
        %{
        if (isempty(qvec))
            break;
        end
        %}
    end
end 

if (~isempty(validname))
    vvec = validoperation(expt.gatherField(validname, varargin{:}));
    qvec = qvec(:,vvec);
end
