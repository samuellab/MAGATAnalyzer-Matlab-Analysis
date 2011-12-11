function magfactor = realUnitsPerPixel (cc)
%function magfactor = realUnitsPerPixel (cc)
%
% how many real units (usually cm) per pixel
% outputs:
%   MAGFACTOR - cm per pixel
% inputs:
%   CC < CameraCalibration

I = randperm(length(cc.camx));
cx = cc.camx(I);
cy = cc.camy(I);
rx = cc.realx(I);
ry = cc.realy(I);

dc = sqrt(diff(cx).^2 + diff(cy).^2);
dr = sqrt(diff(rx).^2 + diff(ry).^2);

mf = dr./dc;
magfactor = median(mf(isfinite(mf)));