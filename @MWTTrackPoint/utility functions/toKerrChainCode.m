function ccstr = toKerrChainCode(dx, dy)
dx = round(dx);
dy = round(dy);
if (any (dx == dy | abs(dx)+abs(dy) > 1))
    error ('need 4 connected contour with single spacing;');
end
mat = [NaN 2 NaN; 0 NaN 1; NaN 3 NaN];

bc = uint8(mat(sub2ind(size(mat),dy+2,dx+2)));
switch (mod(length(dx), 3))
    case 2
        bc = [bc 0];
    case 1
        bc = [bc 0 1];
end
ccstr = char(((bitshift(bc(1:3:end),4) + bitshift(bc(2:3:end),2) + bc(3:3:end))) + '0');
