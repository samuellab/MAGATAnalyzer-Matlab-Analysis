function [pt,inPt1Pt2,inPt3Pt4] = intersectionPoint(pt1,pt2,pt3,pt4)
%function pt = intersectionPoint(pt1,pt2,pt3,pt4)
%
%returns the location of the intersection of the two lines defined by
%pt1-pt2 and pt3-pt4


d = (pt4(2)-pt3(2))*(pt2(1)-pt1(1)) - (pt4(1)-pt3(1))*(pt2(2)-pt1(2));
n = (pt4(1)-pt3(1))*(pt1(2)-pt3(2)) - (pt4(2)-pt3(2))*(pt1(1)-pt3(1));
n2 = (pt2(1)-pt1(1))*(pt1(2)-pt3(2)) - (pt2(2)-pt1(2))*(pt1(1)-pt3(1));
pt = pt1 + n/d*(pt2-pt1);

inPt1Pt2 = 0 <= n/d && n/d <=1;
inPt3Pt4 = 0  <= n2/d && n2/d <=1;