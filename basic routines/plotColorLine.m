function h = plotColorLine (x, y, z, cmap, zrange, varargin)
%function h = plotColorLine (x, y, z, cmap, zrange, varargin)

showcolorbar = false;
varargin = assignApplicable(varargin);

if (~exist('zrange','var') || isempty(zrange))
    zrange = [min(z) max(z)];
end
ih = ishold;
cx = ((1:length(cmap)) - 1)*(zrange(2)-zrange(1))/length(cmap) + zrange(1);
for j = 1:length(x)
    if (isempty(varargin))
        h(j) = plot(x(j),y(j),'b.'); hold on
    else
        h(j) = plot(x(j),y(j),'b.',varargin{:}); hold on
    end
    set(h(j),'Color',interp1(cx,cmap,z(j),'nearest','extrap'));
end

if (showcolorbar)
    ch = colorbar ('vert');%,'CLim',[1 length(cmap)]); 
    yr = get(ch, 'YLim');
    yt = get(ch, 'YTick');

    m = (diff(zrange))/diff(yr);
    b = zrange(1) - m*yr(1);
    zt = m*yt + b;
    for j = 1:length(yt)
        ztl{j} = num2str(zt(j));
    end
    set (ch, 'YTickLabel', ztl);
 
end
if (~ih)
    hold off;
end