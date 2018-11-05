function k = radialDistortion (xp, yp, xd, yd)
%function k = radialDistortion (xp, yp, xd, yd)
%
%xp and yp are the projected points without distortion
%xd and yd are the actual distorted points
%
%we solve for k, s.t. (xp,yp) = (1 + k(1)r^2 + k(2)r^4 + k(3)r^6)*(xd,yd)

x = xp./xd - 1;
y = yp./yd - 1;
r = xp.^2 + yp.^2;

A = [r;r.^2;r.^3];
%size(A)

X = [x y];
%size(X)
AA = [A A];
%size(AA)
%size(pinv(AA))
k = X*pinv(AA);