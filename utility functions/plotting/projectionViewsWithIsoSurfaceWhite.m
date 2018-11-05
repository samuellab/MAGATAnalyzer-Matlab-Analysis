function handles = projectionViewsWithIsoSurfaceWhite(pvstruct, clim)
%function  [ax1, ax2,ax3, ax4] = projectionViews(pvstruct, clim)

clf(gcf);
existsAndDefault('xscale', 1);
existsAndDefault('yscale', 1);
existsAndDefault('zscale', 1);

xaxis = pvstruct.xaxis;
yaxis = pvstruct.yaxis;
zaxis = pvstruct.zaxis;

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

ax(1) = axes('position', [sp1 1-sp2 - ys/h xs/w ys/h]);
ax(2) = axes('position', [1-sp1-zs/w 1-sp2 - ys/h zs/w ys/h]);
ax(3) = axes('position', [sp1 sp2 xs/w zs/h]);
ax(4) = axes('position', [1-sp1-zs/w sp2 zs/w zs/h]);

if (~exist('clim', 'var') || isempty(clim))
    clim = [0 max(pvstruct.imzp(:))];
end

imagesc(xaxis,yaxis,pvstruct.imzp,  'Parent', ax(1));
imagesc(zaxis,yaxis,pvstruct.imxp,  'Parent', ax(2));
imagesc(xaxis,zaxis,pvstruct.imyp,  'Parent', ax(3));

hold(ax(1), 'on');
[~,handles.cont_h(1)] = contour(ax(1), xaxis, yaxis, pvstruct.imzp, [pvstruct.isoval pvstruct.isoval], 'k--', 'LineWidth', 2);
hold(ax(1), 'off');

hold(ax(2), 'on');
[~,handles.cont_h(2)] = contour(ax(2), zaxis, yaxis, pvstruct.imxp, [pvstruct.isoval pvstruct.isoval], 'k--', 'LineWidth', 2);
hold(ax(2), 'off');

hold(ax(3), 'on');
[~,handles.cont_h(3)] = contour(ax(3), xaxis, zaxis, pvstruct.imyp, [pvstruct.isoval pvstruct.isoval], 'k--', 'LineWidth', 2);
hold(ax(3), 'off');


set(ax(1), 'YAxisLocation', 'left', 'XAxisLocation', 'top', 'box', 'on', 'Clim', clim);
set(ax(2), 'YAxisLocation', 'right', 'XAxisLocation', 'top', 'box', 'on', 'Clim', clim,'XDir','reverse');
set(ax(3), 'YAxisLocation', 'left', 'XAxisLocation', 'bottom', 'box', 'on', 'Clim', clim);

axis(ax(1:3), 'equal');
axis(ax(1:3), 'xy');
p = patch(pvstruct.fv, 'parent', ax(4));
set(p,'FaceColor','red','EdgeColor','none','VertexNormals', pvstruct.ison);
daspect(ax(4), [1,1,1])
view(ax(4), 3);
axes(ax(4));

axis(ax(4), pvstruct.axislims);
hlight = camlight;
lighting phong;

handles.ax = ax;
handles.p = p;
handles.light = hlight;

set(gcf, 'color', 'w'); set(ax, 'color', 'w', 'XColor', 'k', 'YColor', 'k', 'ZColor', 'k');
set(ax(4), 'XTick', [], 'YTick', [], 'ZTick', [], 'XColor', 'w', 'YColor', 'w', 'ZColor', 'w'); 
set(ax(2), 'XDir','reverse');
colormap hot;
% 
% axis(ax1, 'equal');
% axis(ax2, 'equal');
% axis(ax3, 'equal');
% axis(ax1, 'xy'); 

