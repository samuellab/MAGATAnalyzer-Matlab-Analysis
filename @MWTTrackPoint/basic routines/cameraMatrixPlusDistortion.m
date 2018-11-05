function [CM, k] = cameraMatrixPlusDistortion(u,v,x,y,z,kstart)
%function [CM, k] = cameraMatrixPlusDistortion(u,v,x,y,z,kstart)
%x,y,z are actual physical locations
%u,v are measured pixel locations relative to the center point of the
%camera
%
%tup,tvp,t = CM*(x,y,z,1)  
%u,v = (up,vp) / (1 + k(1)*r^2 + k(2)*r^4 + k(3)*r^6)

niter = 1000;

if (~exist('kstart','var') || isempty(kstart))
    k = [0 0 0];
else
    k = kstart;
end
for nn = 1:niter
    r = u.^2 + v.^2;
    a = (1 + k(1)*r + k(2)*r.^2 + k(3)*r.^3);
    CM = cameraMatrix(u.*a,v.*a,x,y,z);
    for j = 1:length(x)
        [uu,vv] = projectPt (x(j),y(j),z(j), CM);
        up(j) = uu;
        vp(j) = vv;
    end
    k = radialDistortion (up, vp, u, v);
    
end


end

function [u,v] = projectPt (x,y,z,CM)

z = CM*([x;y;z;1]);
u = z(1)/z(3);
v = z(2)/z(3);

end
