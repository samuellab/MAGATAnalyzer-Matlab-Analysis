function h = filledRibbon (x,y,z,u,v,w,c, varargin)
%function filledRibbon (x,y,z,u,v,w,c, varargin)
%
%plots a ribbon spanning the area between the lines x,y,z and x+u,y+v,z+w
%in the color c
%varargin is passed directly to patch
%returns a handle to the patch graphic created


if ~exist('w', 'var') || isempty(w)
    w = 0;
end
if ~exist('u', 'var') || isempty(u)
    u = 0;
end
if ~exist('v', 'var') || isempty(v)
    v = 0;
end
if ~exist('c', 'var') || isempty(c)
    c = 'b';
end


%make all vectors 1xN 
x = reshape(x,1,[]);
y = reshape(y,1,[]);
z = reshape(z,1,[]);

%if any offsets are scalar, expand to a vector
if all(size(u) == 1)
    u = repmat(u, size(x));
end

if all(size(v) == 1)
    v = repmat(v, size(x));
end
if all(size(w) == 1)
    w = repmat(w, size(x));
end

%make up a set of regions that span the space between the lines

xr = [x(1:end-1); x(1:end-1) + u(1:end-1); x(2:end) + u(2:end); x(2:end)];
yr = [y(1:end-1); y(1:end-1) + v(1:end-1); y(2:end) + v(2:end); y(2:end)];
zr = [z(1:end-1); z(1:end-1) + w(1:end-1); z(2:end) + w(2:end); z(2:end)];

%plot the regions with no edges
h = patch(xr,yr,zr,c, 'LineStyle','none', varargin{:});