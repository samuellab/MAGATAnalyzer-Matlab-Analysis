function [gq, pt0, x, y, theta] = lightDirectionFromImage (rawim, cc, varargin)
% function [gq, pt0, x, y, theta] = lightDirectionFromImage (rawim, cc, varargin)
%
% creates a global quantity that tells the direction to the light source at
% any point;  given a pin shadow image
% if this breaks or gives bad results, try working through each subfunction
% call
%
% inputs: rawim - an image of the pins casting shadows
%         cc < CameraCalibration - camera calibration that converts image to real
%         coords
% outputs: gq - adds field angToLightSource
%          pt0 - the x,y, location of the light source
%          x,y,theta the measured location and angles to light source (at
%          ends of shadows)

predominantTheta = 0;
minLen = 1;
varargin = assignApplicable(varargin);
[im,rx,ry] = cc.morphCamToReal(rawim);
lines = lineEdgesInImage(im);
for j = 1:length(lines)
    lines(j).x = interp1(rx, lines(j).x);
    lines(j).y = interp1(ry, lines(j).y);
    lines(j).len = sqrt(diff(lines(j).x).^2 + diff(lines(j).y).^2);
end

[pt0,x,y,theta] = lightDirectionFromLines(lines, predominantTheta, minLen);

gq = GlobalQuantity();
gq.fieldname = 'angToLightSource';
gq.xField = 'iloc';
gq.xData = pt0;
gq.derivationMethod = @(xin, xData, yData) atan2(xData(2)-xin(2,:), xData(1)-xin(1,:));
