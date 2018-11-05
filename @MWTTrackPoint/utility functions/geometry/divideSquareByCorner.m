function [a,pts] = divideSquareByCorner (ll, ur, x0, x1, x2)
%function [a,pts] = divideSquareByCorner (ll, ur, x0, x1, x2)
%
%given a line between x0 and x1 and a line between x1 and x2, with x1 in the square
%defined by corners at ll, ur; ll < ur
%
%finds the intersection points A and B between the lines x0-x1 (A) and
%x1-x2 (B) and the square
%
%gives the area of the polygon A, (0-3 corners of square), B.  The corners
%are chosen so that the line has a counter-clockwise orientation.

corners = [[ll(1);ll(2)] [ur(1);ll(2)] [ur(1);ur(2)] [ll(1);ur(2)]];

%find points where lines intersect sides
for c1 = 1:4
    k = mod(c1,4)+1;
    [A,in1,in2] = intersectionPoint(corners(:,c1), corners(:,k), x0, x1);
    if (in1 && in2)
        break;
    end
end

for c2 = 1:4
    k = mod(c2,4)+1;
    [B,in1,in2] = intersectionPoint(corners(:,c2), corners(:,k), x1, x2);
    if (in1 && in2)
        break;
    end
end
if (c1 == c2) %intersection through the same side
    pts = [A x1 B];
    threepts = true;
else
    threepts = false;
    %make a set of points oriented counterclockwise
    %A x1 B c2:c1
    dx1 = A - corners(:,c1); dx2 = x1 - A; 
    if (dx1(1)*dx2(2)-dx2(1)*dx1(2)) < 0
        c1 = mod(c1,4)+1;
    end
    
    dx1 = B-x1; dx2 = corners(:,c2) - B; 
    if (dx1(1)*dx2(2)-dx2(1)*dx1(2)) < 0
        c2 = mod(c2,4)+1;
    end
    
    if (c2 <= c1)
        pts = [A x1 B corners(:,c2:c1)];
    else
        pts = [A x1 B corners(:,[c2:4 1:c1])];
    end
end
    %{
plot (corners(1,:), corners(2,:), 'k.-', corners(1,[4 1]), corners(2,[4 1]), 'm-', [x0(1), x1(1), x2(1)], [x0(2), x1(2), x2(2)], 'g.-', pts(1,:), pts(2,:), 'co', corners(1,c1), corners(2,c1), 'r*', corners(1,c2), corners(2,c2), 'g*'); 
text (A(1), A(2), 'A'); text(B(1), B(2), 'B'); text (x0(1), x0(2), 'x0'); text (x1(1), x1(2), 'x1');text (x2(1), x2(2), 'x2'); 
pause
%}
a = polyarea_signed(pts(1,:), pts(2,:));
if (a < 0)
    
    if (~threepts)
        warning ('DSQC: produced points in wrong orientation');
%         plot (corners(1,:), corners(2,:), 'k.-', corners(1,[4 1]), corners(2,[4 1]), 'm-', [x0(1), x1(1), x2(1)], [x0(2), x1(2), x2(2)], 'g.-', pts(1,:), pts(2,:), 'co', corners(1,c1), corners(2,c1), 'r*', corners(1,c2), corners(2,c2), 'g*'); 
%         text (A(1), A(2), 'A'); text(B(1), B(2), 'B'); text (x0(1), x0(2), 'x0'); text (x1(1), x1(2), 'x1');text (x2(1), x2(2), 'x2'); 
%         pause
        a = -a;
    else
        a = (ur(2)-ll(2))*(ur(1)-ll(1)) + a;
    end
end
