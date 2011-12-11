function u = realPtsToCamera(x, camcalinfo)
%function u = realPtsToCamera(x, camcalinfo)
%
%x is a 2 X N vector of points
%we assume z = 0
%camera calibration is given in camcalinfo
%if camcalinfo is [] or nonexistent, u = x
%otherwise, u is the set of camera locations that corresponds to x in real
%space

if (~exist('camcalinfo', 'var') || isempty(camcalinfo))
    u = x;
    return;
end

[ux,uy] = cameraPointFromRealPoint(x(1,:), x(2,:), zeros(size(x(1,:))), camcalinfo);
u = [ux;uy];
