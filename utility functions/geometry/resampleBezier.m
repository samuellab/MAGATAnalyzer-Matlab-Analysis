function q = resampleBezier(q, nsamp)
%function q = resampleBezier(q, nsamp)
%finds best representation of same curve, given restriction that all points
%are equally spaced;  this may just create a line or some crap like that
%
if (nargin < 2)
    nsamp = 100;
end
s = linspace(0,1,nsamp);
M = bezierMatrix(size(q,2)-1, s);

u = q*M;
l = [0 cumsum(sqrt(sum(diff(u,[],2).^2)))];
u = interp1(l/l(end),u',s, 'linear')';

q = (pinv(M')*u')';
