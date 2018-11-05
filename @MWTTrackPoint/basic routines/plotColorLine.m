function hh = plotColorLine (x, y, z, cmap, zrange, varargin)
%function h = plotColorLine (x, y, z, cmap, zrange, varargin)

showcolorbar = false;
decimate = length(x) > 1E4;
decimateLength = 1e4;
varargin = assignApplicable(varargin);


x = interp1(find(isfinite(x)), x(isfinite(x)), 1:length(x), 'nearest', 'extrap');
y = interp1(find(isfinite(y)), y(isfinite(y)), 1:length(y), 'nearest', 'extrap');
z = interp1(find(isfinite(z)), z(isfinite(z)), 1:length(z), 'nearest', 'extrap');


if decimate
    df = round(length(x) / decimateLength);
    k = ones([1 df]);
    d = conv(ones(size(x)), k, 'same');
    x = conv(x, k, 'same')./d;
    y = conv(y, k, 'same')./d;
    z = conv(z, k, 'same')./d;
    x = x(1:df:end);
    y = y(1:df:end);
    z = z(1:df:end);
end
    

existsAndDefault('cmap', jet(256));
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
if (nargout > 0)
    hh = h;
end