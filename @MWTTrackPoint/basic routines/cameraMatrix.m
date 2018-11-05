function CM = cameraMatrix (u,v,x,y,z)
%function CM = cameraMatrix (u,v,x,y,z)
%
%calculates the camera matrix that transforms a point in x,y,z real space
%to a point on the ccd u,v .  x,y,z are a list of real points while u,v are
%the locations of the corresponding points on the CCD sensor.
%
%The camera matrix operates using homogenous coordinates, thus 
%(tu,tv,t) = CM*(x,y,z,1) for some value of t.

%solve P*CM=(u1,v1,u2,v2 . . .)
%where P = [(x1 y1 z1 1 0 0 0 0 -u1x1 -u1y1 -u1z1),(0 0 0 0 x1 y1 z1 1
%-v1x1 -v1y1 -v1z1), . . . ]
%and CM is [c11 c12 c13 c14 c23 . . . c33]
%CM34 = 1, using up 1 free scaling factor

P = zeros(11,length(u)*2);
X = zeros(length(u)*2,1);
for j = 1:length(u)
    k = 2*j - 1;
    P(:,k) = [x(j) y(j) z(j) 1 0 0 0 0 -u(j)*x(j) -u(j)*y(j) -u(j)*z(j)];
    P(:,k+1) = [0 0 0 0 x(j) y(j) z(j) 1 -v(j)*x(j) -v(j)*y(j) -v(j)*z(j)];
    X(k) = u(j);
    X(k+1) = v(j);
end
%size(X)
%size(P')
%size(pinv(P'))
CM = pinv(P')*X;
CM(end+1) = 1;
CM = reshape(CM,4,3);
CM = CM';
    

