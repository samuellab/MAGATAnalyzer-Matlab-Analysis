function [ax1, ax2,ax3] = projectionViews(im3d, xscale, yscale, zscale, clim)
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

%figure out how big to make axes
xf = 1;
yf = diff(yaxis([1 end]))/diff(xaxis([1 end]));
zf = diff(zaxis([1 end]))/diff(xaxis([1 end]));

p = get(gcf, 'position');
w = p(3); h = p(4);

pixr = min(w/(xf + zf + 0.2), h/(yf + zf + 0.2));

xs = pixr*xf;
ys = pixr*yf;
zs = pixr*zf;

sp1 = (w - (xs + zs))/(2.5*w);
sp2 = (h - (ys + zs))/(2.5*h);

ax1 = axes('position', [sp1 1-sp2 - ys/h xs/w ys/h]);
ax2 = axes('position', [1-sp1-zs/w 1-sp2 - ys/h zs/w ys/h]);
ax3 = axes('position', [sp1 sp2 xs/w zs/h]);

if (~exist('clim', 'var') || isempty(clim))
    clim = [0 max(im3d(:))];
end

imagesc(xaxis,yaxis,squeeze(max(im3d,[],3)),  'Parent', ax1);
imagesc(zaxis,yaxis,squeeze(max(im3d,[],2)),  'Parent', ax2);
imagesc(xaxis,zaxis,squeeze(max(im3d,[],1))',  'Parent', ax3);

set(ax1, 'YAxisLocation', 'left', 'XAxisLocation', 'top', 'box', 'on', 'Clim', clim);
set(ax2, 'YAxisLocation', 'right', 'XAxisLocation', 'top', 'box', 'on', 'Clim', clim);
set(ax3, 'YAxisLocation', 'left', 'XAxisLocation', 'bottom', 'box', 'on', 'Clim', clim);


axis(ax1, 'equal');
axis(ax2, 'equal');
axis(ax3, 'equal');

