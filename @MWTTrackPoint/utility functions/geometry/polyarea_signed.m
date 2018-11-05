function area = polyarea_signed(x,y,dim)
%POLYAREA_SIGNED Area of polygon with orientation taken into account.
%   POLYAREA(X,Y) returns the area of the polygon specified by
%   the vertices in the vectors X and Y.  If X and Y are matrices
%   of the same size, then POLYAREA returns the area of
%   polygons defined by the columns X and Y.  If X and Y are
%   arrays, POLYAREA returns the area of the polygons in the
%   first non-singleton dimension of X and Y.  
%
%   The polygon edges must not intersect.  If they do, POLYAREA
%   returns the difference between the counterclockwise
%   encircled areas and the clockwise encircled areas.
%
%   area is signed and > 0 if the points are oriented counterclockwise
%   area is < 0 if the points are oriented clockwise
%
%   POLYAREA(X,Y,DIM) returns the area of the polygons specified
%   by the vertices in the dimension DIM.
%
%   Class support for inputs X,Y:
%      float: double, single

%   Copyright 1984-2004 The MathWorks, Inc.
%   $Revision: 1.12.4.2 $  $Date: 2004/03/02 21:47:55 $

if nargin==1 
  error('MATLAB:polyarea:NotEnoughInputs', 'Not enough inputs.'); 
end

if ~isequal(size(x),size(y)) 
  error('MATLAB:polyarea:XYSizeMismatch', 'X and Y must be the same size.'); 
end

if nargin==2,
  [x,nshifts] = shiftdim(x);
  y = shiftdim(y);
elseif nargin==3,
  perm = [dim:max(length(size(x)),dim) 1:dim-1];
  x = permute(x,perm);
  y = permute(y,perm);
end

siz = size(x);
if ~isempty(x),
  area = -reshape(sum( (x([2:siz(1) 1],:) - x(:,:)).* ...
                 (y([2:siz(1) 1],:) + y(:,:)))/2,[1 siz(2:end)]);
else
  area = sum(x); % SUM produces the right value for all empty cases
end

if nargin==2,
  area = shiftdim(area,-nshifts);
elseif nargin==3,
  area = ipermute(area,perm);
end
