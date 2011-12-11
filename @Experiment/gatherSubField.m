function qvec = gatherSubField (expt, field, subfield, varargin)
%gets all track.(field).(subfield)
%function qvec = gatherSubField (expt, field, subfield, varargin)
%
%outputs:
%QVEC: a kxN array where k is the dimension of field.subfield and N is the
%   number of points
%
%inputs:
%EXPT: a member of the experiment class
%FIELD: a property of EXPT.track
%   if field is 'firsths', gathers the first headsweep from every
%   reorientation only
%   if field is 'lasths', gathers the last headsweep from every reo
%SUBFIELD: a property of EXPT.track.FIELD
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

expandToInds = false;
varargin = assignApplicable(varargin);




if (strcmpi (field, 'firsths'))
    r = [expt.track.reorientation];
    r = r([r.numHS] > 0);
    f = repmat(HeadSwing, size(r));
    for k = 1:length(r)
        f(k) = r(k).headSwing(1);
    end
else
    if (strcmpi (field, 'lasths'))
        r = [expt.track.reorientation];
        r = r([r.numHS] > 0);
        f = repmat(HeadSwing, size(r));
        for k = 1:length(r)
            f(k) = r(k).headSwing(end);
        end
    else
        f = [expt.track.(field)];
    end
end
qvec = [f.(subfield)];

if (expandToInds)
    inds = [f.inds];
    k = 0;
    for j = 1:length(f)
        n = length(f(j).inds);
        inds(k + (1:n)) = j;
        k = k+n;
    end
    qvec = qvec(:,inds);
end