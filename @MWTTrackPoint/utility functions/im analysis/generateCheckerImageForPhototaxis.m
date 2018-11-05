function gradientImage = generateCheckerImageForPhototaxis(x1, x2, y1, y2, cc, checkerSizeInCm, borderSizeInCm, varargin)
%function gradientImage = generateGradientImageForPhototaxis(plateIm, cc, borderSizeInCm, gradstyle, varargin)

existsAndDefault('borderSizeInCm' , 0);

rx = linspace(min(cc.realx), max(cc.realx), 1280);
ry = linspace(min(cc.realy), max(cc.realy), 800);

x0 = (min(x1(1),x2(1)));
w = (max(x2(1),x1(1))) - x0;

xc = mean([x1 x2]);
yc = mean([y1 y2]);
[px,py] = meshgrid(rx,ry);
gi = 255*xor(mod(floor((px-xc) / checkerSizeInCm),2) , mod(floor((py-yc) / checkerSizeInCm),2));
%imagesc(rx,ry,gi); pause
%set region outside [x1,x2,y1,y2] to 0
%set region borderSizeInCm within to 255
gi (px < min([x1 x2]) + borderSizeInCm) = 255;
gi (py < min([y1 y2]) + borderSizeInCm) = 255;
gi (px > max([x1 x2]) - borderSizeInCm) = 255;
gi (py > max([y1 y2]) - borderSizeInCm) = 255;
gi (px < min([x1 x2])) = 0;
gi (py < min([y1 y2])) = 0;
gi (px > max([x1 x2])) = 0;
gi (py > max([y1 y2])) = 0;


gradientImage = cc.morphRealToCam(gi, rx, ry, 'camxaxis', min(cc.camx):max(cc.camx), 'camyaxis', min(cc.camy):max(cc.camy));
gradientImage = uint8(gradientImage);
image(gradientImage); colormap gray(256);
title ('checker');
pause(0.1);