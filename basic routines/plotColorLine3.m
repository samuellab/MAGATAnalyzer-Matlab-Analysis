function hh = plotColorLine3 (x, y, z, u, cmap, urange, varargin)
%function h = plotColorLine (x, y, z, u, cmap, urange, varargin)

showcolorbar = false;
decimate = length(x) > 1E4;
decimateLength = 1e4;
varargin = assignApplicable(varargin);


x = interp1(find(isfinite(x)), x(isfinite(x)), 1:length(x), 'nearest', 'extrap');
y = interp1(find(isfinite(y)), y(isfinite(y)), 1:length(y), 'nearest', 'extrap');
z = interp1(find(isfinite(z)), z(isfinite(z)), 1:length(z), 'nearest', 'extrap');
u = interp1(find(isfinite(u)), u(isfinite(u)), 1:length(u), 'nearest', 'extrap');


if decimate
    df = round(length(x) / decimateLength);
    k = ones([1 df]);
    d = conv(ones(size(x)), k, 'same');
    x = conv(x, k, 'same')./d;
    y = conv(y, k, 'same')./d;
    z = conv(z, k, 'same')./d;
    u = conv(z, k, 'same')./d;
    x = x(1:df:end);
    y = y(1:df:end);
    z = z(1:df:end);
    u = u(1:df:end);
end
    

existsAndDefault('cmap', jet(256));
if (~exist('urange','var') || isempty(urange))
    urange = [min(u) max(u)];
end
ih = ishold;
cx = ((1:length(cmap)) - 1)*(urange(2)-urange(1))/length(cmap) + urange(1);



for j = 1:length(x)
    if (isempty(varargin))
        h(j) = plot3(x(j),y(j),z(j), 'b.'); hold on
    else
        h(j) = plot(x(j),y(j),z(j), 'b.',varargin{:}); hold on
    end
    set(h(j),'Color',interp1(cx,cmap,u(j),'nearest','extrap'));
end

if (showcolorbar)
    ch = colorbar ('vert');%,'CLim',[1 length(cmap)]); 
    yr = get(ch, 'YLim');
    yt = get(ch, 'YTick');

    m = (diff(urange))/diff(yr);
    b = urange(1) - m*yr(1);
    zt = m*yt + b;
    for j = 1:length(yt)
        utl{j} = num2str(ut(j));
    end
    set (ch, 'YTickLabel', utl);
 
end
if (~ih)
    hold off;
end
if (nargout > 0)
    hh = h;
end