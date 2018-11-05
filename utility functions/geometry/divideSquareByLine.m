function [a,pts] = divideSquareByLine (ll, ur, x0, x1)
%function divideSquareByLine (ll, ur, x0, x1)
%
%given a line between x0 and x1 and a square with corners at ll, ur; 
%ll < ur
%
%finds the intersection points A and B between the line and the square
%A is closer to x0 than B
%
%gives the area of the polygon A, (1-3 corners of square), B.  The corners
%are chosen so that the line has a counter-clockwise orientation.

corners = [[ll(1);ll(2)] [ur(1);ll(2)] [ur(1);ur(2)] [ll(1);ur(2)]];

%find corners that are opposite each other

dx = x1 - x0;
pdx = repmat([-dx(2); dx(1)],1,4);

c = corners - [x0(:) x0(:) x0(:) x0(:)];

s = sign(dot(pdx, c));

if (all(sign(s) >= 0))
    a = (ur(2)-ll(2))*(ur(1)-ll(1));
    return;
end
if (all(sign(s) <= 0))
    a = 0;
    return;
end

cp = find(diff(s([1:end 1])) ~= 0);
if (length(cp) ~= 2)
    warning ('did not find exactly two crossing points!');
end

c1 = corners(:,cp(1));
c2 = corners(:,cp(1)+1);

c3 = corners(:,cp(2));
c4 = corners(:,mod(cp(2),4)+1);

pt1 = intersectionPoint(c1, c2, x0, x1);
pt2 = intersectionPoint(c3, c4, x0, x1);

if (dot(dx, pt2-pt1) > 0)
    A = pt1;
    B = pt2;
    c1 = cp(1);
    c2 = cp(2);
else
    A = pt2;
    B = pt1;
    c2 = cp(1);
    c1 = cp(2);
end

 dx1 = A - corners(:,c1); dx2 = B - A; 
 if (dx1(1)*dx2(2)-dx2(1)*dx1(2)) < 0
     c1 = mod(c1,4)+1;
 end
 
 dx1 = B-A; dx2 = corners(:,c2) - B;
 if (dx1(1)*dx2(2)-dx2(1)*dx1(2)) < 0
     c2 = mod(c2,4)+1;
 end
 
 if (c2 <= c1)
     pts = [A B corners(:,c2:c1)];
 else
     pts = [A B corners(:,[c2:4 1:c1])];
 end
    
    
% 
% if (dot(dx, pt2-pt1) > 0)
%     pts = [pt1(:) corners(:,s>0) pt2(:)];
% else
%     pts = [pt2(:) corners(:,s>0) pt1(:)];
% end
%pcolor(im); hold on; axis equal;
a = polyarea(pts(1,:), pts(2,:));
%{
plot (corners(1,:), corners(2,:), 'k.-', corners(1,[4 1]), corners(2,[4 1]), 'm-', [x0(1), x1(1)], [x0(2), x1(2)], 'g-',x0(1),x0(2),'g*', pts(1,:), pts(2,:), 'co'); 
hold on;
patch(pts(1,:), pts(2,:), a);
text (A(1), A(2), 'A'); text(B(1), B(2), 'B'); text (x0(1), x0(2), 'x0'); text (x1(1), x1(2), 'x1');
%}
%hold off;