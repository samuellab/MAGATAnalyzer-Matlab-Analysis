function pvstruct = projectionViewsWithIsoSurfaceCalculations(im3d, xscale, yscale, zscale, isoval)
%function projectionViews(im3d, xaxis, yaxis, zaxis, clim)
%function projectionViews(im3d, xscale, yscale, zscale, clim)

clf(gcf);
existsAndDefault('xscale', 1);
existsAndDefault('yscale', 1);
existsAndDefault('zscale', 1);

if (length(xscale) == 1)
    xaxis = (1:size(im3d,2))*xscale;
else
    xaxis = xscale;
end
if (length(yscale) == 1)
    yaxis = (1:size(im3d,1))*yscale;
else
    yaxis = yscale;
end
if (length(zscale) == 1)
    zaxis = (1:size(im3d,3))*zscale;
else
    zaxis = zscale;
end

pvstruct.xaxis = xaxis;
pvstruct.yaxis = yaxis;
pvstruct.zaxis = zaxis;
pvstruct.imzp = squeeze(max(im3d,[],3));
pvstruct.imxp = squeeze(max(im3d,[],2));
pvstruct.imyp = squeeze(max(im3d,[],1))';
pvstruct.fv = isosurface(xaxis, yaxis,zaxis,im3d, isoval);
pvstruct.ison = isonormals(xaxis, yaxis,zaxis,im3d,pvstruct.fv.vertices);
pvstruct.im3d = im3d;
pvstruct.isoval = isoval;

xinds = [find(any(pvstruct.imzp >= isoval, 1), 1, 'first') find(any(pvstruct.imzp >= isoval, 1), 1, 'last')];
yinds = [find(any(pvstruct.imzp >= isoval, 2), 1, 'first') find(any(pvstruct.imzp >= isoval, 2), 1, 'last')];
zinds = [find(any(pvstruct.imxp >= isoval, 1), 1, 'first') find(any(pvstruct.imxp >= isoval, 1), 1, 'last')];

pvstruct.axislims = [xaxis(xinds) yaxis(yinds) zaxis(zinds)];