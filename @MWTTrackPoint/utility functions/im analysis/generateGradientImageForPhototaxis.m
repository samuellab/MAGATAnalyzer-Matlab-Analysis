function gradientImage = generateGradientImageForPhototaxis(x1, x2, cc, borderSizeInCm, gradstyle, varargin)
%function gradientImage = generateGradientImageForPhototaxis(plateIm, cc, borderSizeInCm, gradstyle, varargin)

existsAndDefault('borderSizeInCm' , 0);
existsAndDefault('gradstyle', 'linear');

rx = linspace(min(cc.realx), max(cc.realx), 1280);
ry = linspace(min(cc.realy), max(cc.realy), 800);

x0 = (min(x1(1),x2(1)));
w = (max(x2(1),x1(1))) - x0;
switch lower(gradstyle)
    case {'linear', 'lin'}
        gx = 255*(rx - x0)/w; 
        %gx(gx > 255) = 0;
    case {'exponential', 'exp'}
        gx = 255.^((rx - x0)/w);
        %gx(gx > 255) = 0;
    case {'power', 'pow'}
        p = varargin{1};
        gx = 255*(((rx - x0)/w).^p);
        %gx(gx > 255) = 0;
    otherwise
        disp('allowed styles are linear, exponenential, power ');
        gradientImage = [];
        return;
end
gi = repmat(gx, [length(ry) 1]);
%inds = find(ry-ry(1) < borderSizeInCm);
gi((ry-ry(1) < borderSizeInCm),:) = 255;
%inds = find(ry(end) - ry < borderSizeInCm);
gi((ry(end) - ry < borderSizeInCm), :) = 255;

gradientImage = cc.morphRealToCam(gi, rx, ry, 'camxaxis', min(cc.camx):max(cc.camx), 'camyaxis', min(cc.camy):max(cc.camy));
%gradientImage = morphImage2(gi, rx, ry, cc.camx, cc.camy, cc.realx, cc.realy, min(cc.camx):max(cc.camx), min(cc.camy):max(cc.camy));
gradientImage(gradientImage > 255) = 0;
gradientImage(gradientImage < 0) = 0;
gradientImage = uint8(gradientImage);
image(gradientImage); colormap gray(256);
title (gradstyle);
pause(0.1);