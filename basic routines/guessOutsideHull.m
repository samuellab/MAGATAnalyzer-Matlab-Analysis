function [srcx, srcy, dstx, dsty] = guessOutsideHull (srcx, srcy, dstx, dsty, dstxaxis, dstyaxis)
%function [srcx, srcy, dstx, dsty] = guessOutsideHull (srcx, srcy, dstx, dsty, dstxaxis, dstyaxis)
if (size(srcx, 1) == 1)
    srcx = srcx';
end
if (size(srcy, 1) == 1)
    srcy = srcy';
end
if (size(dstx, 1) == 1)
    dstx = dstx';
end
if (size(dsty, 1) == 1)
    dsty = dsty';
end

[dx,dy] = avggridspacing(dstx', dsty');
dstxaxis = linspace(min(dstxaxis),max(dstxaxis), ceil((max(dstxaxis) - min(dstxaxis))/dx));
dstyaxis = linspace(min(dstyaxis),max(dstyaxis), ceil((max(dstyaxis) - min(dstyaxis))/dy));
[x,y] = meshgrid(dstxaxis, dstyaxis);
x = x(:);
y = y(:);
npts = ceil(max((max(dstx) - min(dstx))/dx, (max(dsty) - min(dsty))/dy)) + 3;

chull = convhull(dstx, dsty);
inp = inpolygon(x,y,dstx(chull), dsty(chull));
x = x(~inp);
y = y(~inp);
if (isempty(x))
    %no points outside convex hull
    return;
end

new_srcx = zeros(size(x));
new_srcy = new_srcx;
for j = 1:length(x)
    [xd, yd, I] = nearestNPoints (x(j), y(j), dstx, dsty, npts);
    [xg,yg] = affineGuess (x(j), y(j), xd, yd, srcx(I), srcy(I));
    new_srcx(j) = xg;
    new_srcy(j) = yg;
   %{ 
    if (any(~isfinite(xd) | ~isfinite(yd) | ~isfinite(xg) | ~isfinite(yg)))
        xd
        yd
        xg
        yg
        x(j)
        y(j)
        npts
    end
    %}
end
srcx = [srcx;new_srcx];
srcy = [srcy;new_srcy];
dstx = [dstx;x];
dsty = [dsty;y];

function [x2g, y2g] = affineGuess (x0, y0, x1, y1, x2, y2)
%function [x2g, y2g] = affineGuess (x0, y0, x1, y1, x2, y2)
%
%find the affine transformation between x1,y1 and x2,y2, then apply it to
%x0,y0
if (size(x1,1) > 1)
    x1 = x1';
end
if (size(x2,1) > 1)
    x2 = x2';
end
if (size(y1,1) > 1)
    y1 = y1';
end
if (size(y2,1) > 1)
    y2 = y2';
end

R = cp2tform([x1;y1]',[x2;y2]', 'projective');

[x2g, y2g] = tformfwd(R, x0, y0);
% 
% x = [x1; y1; ones(size(x1))];
% y = [x2; y2; ones(size(x1))];
% %y = A x
% %A = y * pinv(x)
% A = y * pinv(x);
% z = A * [x0;y0;1];
% x2g = z(1)/(z(3)+eps);
% y2g = z(2)/(z(3)+eps);
%     
