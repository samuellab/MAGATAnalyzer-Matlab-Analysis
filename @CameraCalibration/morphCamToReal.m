function [realim,realxaxis,realyaxis] = morphCamToReal(cc, camim, varargin)
% function [realim,realxaxis,realyaxis] = morphCamToReal(cc, camim, varargin)
% morphs an image taken by the camera to real space -- 
% pcolor (realxaxis, realyaxis, realim) will produce an image as it appears
%       in real space with appropriate measurements on axes
% outputs:
%   REALIM - a MxN image
%   REALXAXIS - a 1xN list of x-locations; if realxaxis is not passed in,
%       realxaxis is determined automatically to span the image and have
%       the same approximate size as the camera image
%   REALYAXIS - a 1xM list of y-locations; if realyaxis is not passed in,
%       realyaxis is determined automatically to span the image and have
%       the same approximate size as the camera image
% inputs:
%   CC < CameraCalibration
%   CAMIM : a HxW image
%
%  optional args: 
%       camxaxis (default 1:W)
%       camyaxis (default 1:H)
%       realxaxis
%       realyaxis
%       realresolution (if passed, realxaxis, realyaxis have spacing of
%           realresolution)
%
camxaxis = 1:size(camim,2);
camyaxis = 1:size(camim,1);

realxaxis = [];
realyaxis = [];
realresolution = [];
varargin = assignApplicable(varargin);
if (isempty(realxaxis) || isempty(realyaxis))
    
    camleft = [repmat(min(camxaxis), size(camyaxis)); camyaxis];
    camright = [repmat(max(camxaxis), size(camyaxis)); camyaxis(end:-1:1)];
    camtop = [camxaxis; repmat(max(camyaxis), size(camxaxis))];
    cambottom = [camxaxis(end:-1:1); repmat(min(camyaxis), size(camxaxis))];
        
    realleft = cc.realPtsFromCamPts(camleft);
    realtop = cc.realPtsFromCamPts(camtop);
    realright = cc.realPtsFromCamPts(camright);
    realbottom = cc.realPtsFromCamPts(cambottom);
    
    box = [realleft, realtop, realright, realbottom];
    geom = polygeom(box(1,:), box(2,:));
    xc = geom(2);
    yc = geom(3);
    
    [~,I] = min(abs(realleft(1,:) - xc));
    x1 = realleft(1,I);
    [~,I] = min(abs(realleft(2,:) - yc));
    y1 = realleft(2,I);
    
    [~,I] = min(abs(realright(1,:) - xc));
    x2 = realright(1,I);
    [~,I] = min(abs(realright(2,:) - yc));
    y2 = realright(2,I);
    
    [~,I] = min(abs(realtop(1,:) - xc));
    x3 = realtop(1,I);
    [~,I] = min(abs(realtop(2,:) - yc));
    y3 = realtop(2,I);
    
    [~,I] = min(abs(realbottom(1,:) - xc));
    x4 = realbottom(1,I);
    [~,I] = min(abs(realbottom(2,:) - yc));
    y4 = realbottom(2,I);
    
    x0 = min([x1 x2 x3 x4]);
    y0 = min([y1 y2 y3 y4]);
    x1 = max([x1 x2 x3 x4]);
    y1 = max([y1 y2 y3 y4]);
    
    
    %{
    camcorners = [min(camxaxis) max(camxaxis) max(camxaxis) min(camxaxis); min(camyaxis) min(camyaxis) max(camyaxis) max(camyaxis)];
    rcorners = cc.realPtsFromCamPts(camcorners);
    x0 = max(min(rcorners(1,:)), min(cc.realx));
    x1 = min(max(rcorners(1,:)), max(cc.realx));
    y0 = max(min(rcorners(2,:)), min(cc.realy));
    y1 = min(max(rcorners(2,:)), max(cc.realy));
    %}
    
    
    if (isempty(realxaxis))
        if (isempty(realresolution))
            realxaxis = linspace(x0, x1, size(camim,2));
        else
            realxaxis = x0:realresolution:x1;
        end
    end
    if (isempty(realyaxis))
         if (isempty(realresolution))
            realyaxis = linspace(y0, y1, size(camim,1));
        else
            realyaxis = y0:realresolution:y1;
        end
    end
end
    
[rxp,ryp] = meshgrid(realxaxis,realyaxis);
rxp = rxp(:);
ryp = ryp(:);

xc = cc.r2cX(rxp, ryp);
yc = cc.r2cY(rxp, ryp);
sum(~isfinite(xc));
sum(~isfinite(yc));

im2 = interp2(camxaxis, camyaxis, double(camim), xc, yc, '*linear');

realim = reshape(im2, [length(realyaxis) length(realxaxis)]);
