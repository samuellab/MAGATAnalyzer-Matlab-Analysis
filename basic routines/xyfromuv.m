function [x,y] = xyfromuv (u, v, C)
%function [x,y] = xyfromuv (u, v, C)

%gets the x,y coordinate (assumes z = 0) in real space for a position x,y
%on the sensor.  Camera matrix is given in C

M = [C(1,1)-C(3,1)*u, C(1,2)-C(3,2)*u; C(2,1)-C(3,1)*v, C(2,2) - C(3,2)*v];
Y = [-C(1,4)+C(3,4)*u;-C(2,4)+C(3,4)*v];
X = inv(M)*Y;
x=X(1);
y=X(2);
