function camim = morphRealToCam(cc, realim, realxaxis, realyaxis, varargin)
%function camin = morphRealToCam(cc, realim, realxaxis, realyaxis, varargin)
% morphs an image defined in real space to appear as it would if the camera had taken
%   a picture of it
%
% outputs:   
%   CAMIM : a HxW image
% inputs:
%   CC < CameraCalibration
%   REALIM - a MxN image
%   REALXAXIS - a 1xN list of x-locations; if realxaxis is not passed in,
%       realxaxis is determined automatically to span the image and have
%       the same approximate size as the camera image
%   REALYAXIS - a 1xM list of y-locations; if realyaxis is not passed in,
%       realyaxis is determined automatically to span the image and have
%       the same approximate size as the camera image
%
%  optional args: 
%       camxaxis - (1xW) default is to span the calibrated range
%       camyaxis - (1xH) default is to span the calibrated range
camxaxis = min(cc.camx):max(cc.camx);
camyaxis = min(cc.camy):max(cc.camy);
%camxaxis = 1:size(camim,2);
%camyaxis = 1:size(camim,1);

varargin = assignApplicable(varargin);

    
[cxp,cyp] = meshgrid(camxaxis,camyaxis);
cxp = cxp(:);
cyp = cyp(:);

xr = cc.c2rX(cxp, cyp);
yr = cc.c2rY(cxp, cyp);

im2 = interp2(realxaxis, realyaxis, double(realim), xr, yr, '*linear');

camim = reshape(im2, [length(camyaxis) length(camxaxis)]);
