function campts = camPtsFromRealPts(cc, realpts)
% function campts = camPtsFromRealPts(cc, realpts)
% map real locations to a points on the camera sensor
% outputs:
%   CAMPTS - a 2xN list of pixel locations
% inputs:
%   CC < CameraCalibration
%   REALPTS - a 2xN list of real points;
if (isempty(realpts))
    campts = [];
    return;
end
if (~isreal(realpts))
    realpts = [real(realpts);imag(realpts)];
end
realpts = double(realpts);
x = cc.r2cX(realpts(1,:), realpts(2,:));
y = cc.r2cY(realpts(1,:), realpts(2,:));
campts = [x;y];