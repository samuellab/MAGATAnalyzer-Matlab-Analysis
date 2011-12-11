function [x,y] = realPointFromCameraPoint (u, v, CM, K)
%function [x,y] = realPointFromCameraPoint (u, v, CM, K)
%
%finds the point (x,y,z=0) that corresponds to the camera point u, v given
%a camera projection matrix CM and distortion parameters in K

A = inv(CM(:,[1 2 4]));
r = u.^2 + v.^2;
a = (1 + K(1)*r + K(2)*r.^2 + K(3) * r.^3);
%{
size(u)
size(v)
size(a)
size(A)
size([a.*u;a.*v;ones(size(u))])
%}
z = A*[a.*u;a.*v;ones(size(u))];
x = z(1,:)./z(3,:);
y = z(2,:)./z(3,:);