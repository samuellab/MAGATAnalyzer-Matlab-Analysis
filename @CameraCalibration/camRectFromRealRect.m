function [outerrect, innerrect, quadrilateral] = camRectFromRealRect(cc, realrect)
%function [outerrect, innerrect, quadrilateral] = camRectFromRealRect(cc, realrect)
%cc < CAMCALINFO
%realrect = [left bottom right top] (x1,y1, x2,y2)
%
%the rectangle defined by real coords is translated into a quadrilateral in camera space
%outterrect is the rectangle that contains this quadrilateral
%innerrect is the rectangle contained within this quadrilateral
%outterrect, innerrect are in the form [left right bottom top]
%quadrilateral is a 2x4 list of corner locations

realrect = realrect(:)';
quadrilateral = cc.camPtsFromRealPts([realrect([1 3 3 1]);realrect([2 2 4 4])]);
x = sort(quadrilateral(1,:)); 
y = sort(quadrilateral(2,:));

outerrect = [x(1) y(1) x(4) y(4)];
innerrect = [x(2) y(2) x(3) y(3)];