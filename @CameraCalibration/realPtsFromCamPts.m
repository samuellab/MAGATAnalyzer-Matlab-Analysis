function realpts = realPtsFromCamPts(cc, campts)
%function realpts = realPtsFromCamPts(cc, campts)
% map real locations from points on the camera sensor
% outputs:
%   REALPTS - a 2xN list of real points;
% inputs:
%   CC < CameraCalibration
%   CAMPTS - a 2xN list of pixel locations
x = cc.c2rX(campts(1,:), campts(2,:));
y = cc.c2rY(campts(1,:), campts(2,:));
realpts = [x;y];