function [x,y] = straightLineEdgeThruPoint (pt, im, xd, yd)
%function [x,y] = straightLineEdgeThruPoint (pt, im, xd, yd)

debug = false;
thresh = 0.15;

existsAndDefault('xd', []);
existsAndDefault('yd', []);
if (isempty(xd) || isempty(yd))
    [xd,yd] = imgradient(im,3);
end
e = xd.^2 + yd.^2;
e = e/interp2(e, pt(1), pt(2), '*linear');
th = atan2(interp2(xd, pt(1), pt(2), '*linear'),interp2(yd, pt(1), pt(2), '*linear'));
left = pt;
right = left;
keepgoing = true;
stepsize = 10;
while (keepgoing)
    keepgoing = false;
    perpx = sin(th)*linspace(-2,2,5);
    perpy = cos(th)*linspace(-2,2,5);
    parx = cos(th)*stepsize;
    pary = sin(th)*stepsize;
    rval = interp2(e, right(1) + parx + perpx, right(2) + pary + perpy, '*linear', -1);
    [rmv, I] = max(rval);
    if (rmv > thresh)
        keepgoing = true;
        right(1) = right(1) + parx + perpx(I);
        right(2) = right(2) + pary + perpy(I);
    end
    lval = interp2(e, left(1) - parx + perpx, left(2) - pary + perpy, '*linear', -1);
    [lmv, I] = max(lval);
    if (lmv > thresh)
        keepgoing = true;
        left(1) = left(1) - parx + perpx(I);
        left(2) = left(2) - pary + perpy(I);
    end
    th = atan2(right(2)-left(2), right(1)-left(1));
    if (debug)
        imagesc(e); axis image; hold on
        plot ([left(1) right(1)], [left(2) right(2)], 'r.-', 'LineWidth', 2, 'MarkerSize', 10); hold off
        pause(0.01);
    end
    l = sqrt((left(1)-right(1))^2 + (left(2)-right(2))^2);
    stepsize = max(stepsize, l/4);
    %{
    maxv = min(1, max(interp2(e, linspace(left(1),right(1), 20), linspace(left(2),right(2),20))))
    e = e/maxv;
    e = e.*(e > thresh);
    %}
end

if (all(left == right))
    x = [];
    y = [];
    return
end

[par,perp] = meshgrid(linspace(-0.1,1.1,100), linspace(-10,10,20));
xx = par(:)*(right(1)-left(1)) + left(1) + perp(:)*sin(th);
yy = par(:)*(right(2)-left(2)) + left(2) + perp(:)*cos(th);
valid = xx >= 1 & xx <= size(e,2) & yy >=1 & yy <= size(e,1);
xx = xx(valid);
yy = yy(valid);
maxv = max(interp2(e, xx, yy, '*linear'));
e = e/maxv - thresh;
miny = floor(min(yy));
minx = floor(min(xx));
mye = e(miny:ceil(max(yy)),minx:ceil(max(xx)));

op = optimset('fminsearch');
%op.maxIter = 20;
fun = @(x) -sumOfE(mye, x);
x = [left(1)-minx+1 right(1)-minx+1 left(2)-miny+1 right(2)-miny+1];
%imagesc(mye); axis image; hold on
%plot (x(1:2), x(3:4), 'r.-', 'LineWidth', 2, 'MarkerSize', 10); hold off

x = fminsearch(fun,x,op);
%hold on; plot (x(1:2), x(3:4), 'g.-', 'LineWidth', 2, 'MarkerSize', 10); hold off
y = [x(3) x(4)] + miny - 1;
x = [x(1) x(2)] + minx - 1;

function y = sumOfE(e, x)
l = ceil(sqrt((x(2)-x(1)).^2 + (x(4)-x(3)).^2));
y = sum(interp2(e, linspace(x(1), x(2), l), linspace(x(3), x(4), l), '*Linear', -1000));

