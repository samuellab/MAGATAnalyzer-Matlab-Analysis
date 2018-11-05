function [xc, r] = circleFromThreePts (x1,x2,x3)
%function [xc, r] = circleFromThreePts (x1,x2,x3)
%
%finds the center and radius of a circle going through the points x1, x2,
%and x3
x1 = reshape(x1,[1 2]);
x2 = reshape(x2,[1 2]);
x3 = reshape(x3,[1 2]);

M = [x1(1) x1(2) 1; x2(1) x2(2) 1; x3(1) x3(2) 1];
pf = -1/det(M);

v = [x1*x1', x2*x2', x3*x3'];

N = M;
N(:,1) = v';
D = pf * det(N);

N = M;
N(:,2) = v';
E = pf * det(N);

N = M;
N(:,3) = v';
F = pf * det(N);
xc = -[D,E]/2;
r = sqrt(-F + D^2/4 + E^2/4);