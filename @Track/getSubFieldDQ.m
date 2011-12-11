function qvec = getSubFieldDQ (track, subfield, quantityName, varargin)
% gets derived quantity from a TrackPart
% function qvec = getSubFieldDQ (track, subfield, quantityName, varargin)
%
% track.subfield must be a TrackPart 
% calls track.subfield(inds).getDerivedQuantity(quantityName,varargin)
%
% outputs:
%   QVEC: a kxN matrix of values, where qvec(:,j) corresponds to the jth
%   point
% inputs: 
%   TRACK: a member of the track class
%   SUBFIELD: the subfield; e.g. 'run', 'reorientation'
%   QUANTITYNAME: the name of the quantity to get
%   VARARGIN:
%   passing 'inds', inds will use a subset of subfield
%   passing 'indsExpression', expression will use expression to generate inds
%
%   note these are indices of the SUBFIELD (e.g. to get from first second and
%       fourth run, subfield = 'run', inds = [1 2 4])
%
%   if indsexpression is 'notlast' or 'notfirst', we gather all but the last
%       or first in the track
%
% examples
% track.getSubFieldDQ('run','speed','position','mean') will get the mean
% speed for all tracks
% track.getSubFieldDQ('reorientation', 'eti', 'indsExpression', ...
% '[track.reorientation.numHS] == 1', 'position', 'start') 
% will get the start time of all reorientations with exactly 1 headsweep
%
% if both inds and indsExpression are passed, we select the intersection
% (but don't do this, come on already)
%
% iff subfield is 'firsths', we take the first headsweep in each
% reorienation only, and we ignore inds, indsExpression
%
% iff subfield is 'lasths', we take the first headsweep in each
% reorienation only, and we ignore inds, indsExpression
%

% 


inds = [];
indsExpression = [];
varargin = assignApplicable(varargin);

if (strcmpi (subfield, 'firsths'))
     r = [track.reorientation];
     r = r([r.numHS] > 0);
     f = repmat(HeadSwing, size(r));
     for k = 1:length(r)
         f(k) = r(k).headSwing(1);
     end
else
    if (strcmpi (subfield, 'lasths'))
        r = [track.reorientation];
        r = r([r.numHS] > 0);
        f = repmat(HeadSwing, size(r));
        for k = 1:length(r)
            f(k) = r(k).headSwing(end);
        end
    else
        f = [track.(subfield)];
        if (isempty(inds))
            inds = 1:length(f);
        end
        if (~isempty(indsExpression))
            switch(lower(indsExpression))
                case 'notfirst'
                    if (length(inds) > 1)
                        inds = inds(2:end);
                    else
                        inds = [];
                    end
                case 'notlast'
                    if (length(inds) > 1)
                        inds = inds(1:(end-1));
                    else
                        inds = [];
                    end
                otherwise
                    goodinds = find(eval(indsExpression));
                    inds = intersect(inds, goodinds);
            end
        end
        f = f(inds);
    end
end
    
if (isempty(f))
    qvec = [];
    return;
end

if (isempty(varargin))
    %if we don't need any additional computation or selection, then this is
    %fine
    qvec = track.getDerivedQuantity(quantityName, false, [f.inds]);
else
    qvec = [];
    for j = 1:length(f)
        qvtemp = f(j).getDerivedQuantity(quantityName, varargin{:});
        qvec = [qvec qvtemp];            
    end
end