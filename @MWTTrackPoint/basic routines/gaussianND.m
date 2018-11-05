function g = gaussianND (x, x0, c)
%function g = gaussianND (x, x0, c)
%
%x is a NDxMpts vector
%x0 is a NDx1 center 
%c = is a NDxND covariance matrix

%dx = diff(x,1,2);
%area = prod(dx);
%area(end+1) = area(end);
%area(1)

x = x - repmat(x0, [1, size(x,2)]);
max(x,[],2)
min(x,[],2)
min(abs(x),[],2)
nd = size(c,1);
g = (1 / sqrt(2*pi*det(c)))^(nd) * exp(-0.5 * dot (x, c\x));
max(g)
%g = g.*area;
