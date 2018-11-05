function [u,v] = cameraPointFromRealPoint (x,y,z,camcalinfo)
%function [u,v] = cameraPointFromRealPoint (x,y,z,camcalinfo)
%


%first calculate the new positions without taking into account barrel
%distortion
%homogenous coords:
%t*[u v 1] = CM * [x y z 1]
x = reshape(x,1,[]);
y = reshape(y,1,[]);
if (isempty(z))
    z = zeros(size(x));
end
t = ones(size(x));
temp = camcalinfo.CM * [x;y;z;t];
u = temp(1,:)./temp(3,:);
v = temp(2,:)./temp(3,:);

%now calculate barrel distortion
%[u',v'] = (1 + k(1)r^2 + k(2)r^4 + k(3)r^6) * [u,v];
r = u.^2 + v.^2;
m = 1./(1 + camcalinfo.K(1)*r + camcalinfo.K(2)*r.^2 + camcalinfo.K(3)*r.^3);

%finally, offset image by center point
u = m.*u + camcalinfo.xc;
v = m.*v + camcalinfo.yc;


%camera -> real
%first offset, then barrel, then inv(CM)

%
%function [x,y] = realPointFromCameraPoint (u, v, CM, K)
%function [x,y] = realPointFromCameraPoint (u, v, CM, K)
%
%finds the point (x,y,z=0) that corresponds to the camera point u, v given
%a camera projection matrix CM and distortion parameters in K

%A = inv(CM(:,[1 2 4]));
%r = u.^2 + v.^2;
%a = (1 + K(1)*r + K(2)*r.^2 + K(3) * r.^3);

%z = A*[a.*u;a.*v;ones(size(u))];
%x = z(1,:)./z(3,:);
%y = z(2,:)./z(3,:);