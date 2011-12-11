function [dirim,distim, interior] = generateBoundaryImage(imsize, structvar)
%generates images of size imsize with the boundaries specified in
%structvar
%
%dirim gives the angle (in radians -pi to pi) of the normal to the boundary
%distim gives the oriented distance to the nearest boundary
%distim < 0 if the point is inside the squares
%interior is a logical array that is true for points inside the rectangle
%
%points that are not on any boundary have a distance of infinity and a 
%direction of -100
%
%about structvar
%the first rule of structvar is you do not talk about structvar
%the second rule of structvar is you do not talk about structvar
%the third rule of structvar is that it is a Nx2 array of points,
%specifying rectangles
%the fourth rule of structvar is the lower left (in image coordinates; upper left in xy coords)
%corner comes first, then the other 3 points are specified in
%counterclockwise (in image coords; clockwise in xy coords) order
%the fifth rule of structvar is you do not talk about structvar
%

dirim = zeros(imsize) - 100;
distim = Inf(imsize);

pt1 = structvar(1:4:end,:);
pt2 = structvar(2:4:end,:);
pt3 = structvar(3:4:end,:);
pt4 = structvar(4:4:end,:);

v = pt2 - pt1;
dir = atan2(-v(:,1),v(:,2));
[distim,dirim] = drawLineDistanceOnImage(distim, pt1, pt2, dirim, dir);
    
v = pt3 - pt2;
dir = atan2(-v(:,1),v(:,2));
[distim,dirim] = drawLineDistanceOnImage(distim, pt2, pt3, dirim, dir);
    
v = pt4 - pt3;
dir = atan2(-v(:,1),v(:,2));
[distim,dirim] = drawLineDistanceOnImage(distim, pt3, pt4, dirim, dir);
    
v = pt1 - pt4;
dir = atan2(-v(:,1),v(:,2));
[distim,dirim] = drawLineDistanceOnImage(distim, pt4, pt1, dirim, dir);
    
interior = false(imsize);
[x,y] = meshgrid(1:size(dirim,2),1:size(dirim,1));
x = reshape(x,1,[]);
y = reshape(y,1,[]);
for j = 0:(length(pt1)-1)
    r = structvar((1:4)+4*j, :);
    interior(inpolygon(x,y,r(:,1), r(:,2))) = true;
end
distim(interior) = -distim(interior);
