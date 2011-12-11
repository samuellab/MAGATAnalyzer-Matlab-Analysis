function [slr,pl] = straightLineRatio (pts, dn)
%function slr = straightLineRatio (pts, dn)
%slr(j) = distance between pts(j+dn) and pts(j) divided by the path length
%distance between them

sz = size(pts);
if (sz(1) > sz(2))
    pts = pts';
end

dl = sqrt(sum(diff(pts,1,2).^2,1));
pl = conv2(dl, ones([1 dn]));
dx = sqrt(sum((pts(:,(1+dn):end)-pts(:,1:(end-dn))).^2,1));
%{
size(dl)
size(pl)
size(dx)
size(pl(dn:(end-dn)))
%}
slr = dx./(pl(dn:(end-dn+1)));
