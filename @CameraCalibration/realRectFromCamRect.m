function [outerrect, innerrect, quadrilateral] = realRectFromCamRect(cc, camrect)
%function [outerrect, innerrect, quadrilateral] = realRectFromCamRect(cc, camrect)
%cc < CAMCALINFO
%camrect = [left bottom right top] (x1,y1, x2,y2)
%
%the rectangle defined by camrect is translated into a quadrilateral in real space
%outterrect is the rectangle that contains this quadrilateral
%innerrect is the rectangle contained within this quadrilateral
%outterrect, innerrect are in the form [left right bottom top]
%quadrilateral is a 2x4 list of corner locations

camrect = camrect(:)';
quadrilateral = cc.realPtsFromCamPts([camrect([1 3 3 1]);camrect([2 2 4 4])]);
x = sort(quadrilateral(1,:)); 
y = sort(quadrilateral(2,:));

outerrect = [x(1) y(1) x(4) y(4)];
innerrect = [x(2) y(2) x(3) y(3)];