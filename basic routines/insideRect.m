function tf = insideRect (rect, ptlist)
%rect is a 4 element vector with [x0 x1 y0 y1]
%ptlist is a 2XN list of points
%tf(j) = x0 <= x(1,j) <= x1 && y0 < x(2,j) < y1

tf = (rect(1) <= ptlist(1,:)) & (rect(2) >= ptlist(1,:)) &...
     (rect(3) <= ptlist(2,:)) & (rect(4) >= ptlist(2,:));
 