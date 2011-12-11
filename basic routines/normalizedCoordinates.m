function [nx,ny] = normalizedCoordinates(x,y,ax)

if ~existsAndDefault('ax', [])
    ax = gca;
end


xl = get(ax, 'XLim');
yl = get(ax, 'YLim');

units = get(ax, 'Units');
set(ax,'Units', 'normalized');
pos = get(ax, 'Position');
set(ax, 'Units', units);
nx = (x-xl(1))/(xl(2)-xl(1)) * pos(3) + pos(1);
ny = (y-yl(1))/(yl(2)-yl(1)) * pos(4) + pos(2);


