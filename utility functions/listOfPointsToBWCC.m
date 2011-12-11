function [cc,d] = listOfPointsToBWCC (points, connectivity)
%function cc = listOfPointsToBWCC (points, connectivity)

existsAndDefault('connectivity', 8);
dx = diff(points(1,:));
dy = diff(points(2,:));
dx = min(abs(dx(dx ~= 0)));
dy = min(abs(dy(dy ~= 0)));
d = max(dx,dy);
d = max(d,0.01);
x = points(1,[1:end 1])/d;
y = points(2,[1:end 1])/d;


x = round(x - min(x)) + 1;
y = round(y - min(y)) + 1;
%xx = [];
%yy = [];
%{
for j = 1:(length(x) - 1)
    np = max(abs((x(j+1) - x(j))), abs((y(j+1)-y(j)))) + 1;
    xx = [xx, round(linspace(x(j),x(j+1), np))];
    yy = [yy, round(linspace(y(j),y(j+1), np))];
end

plot (x,y,'r.', xx, yy)
%}        
sz = [max(y), max(x)];
[xx,yy] = meshgrid(1:sz(2), 1:sz(1));
[~,on] = inpolygon(xx,yy,x,y);
ind = sub2ind(sz, yy(on), xx(on));
cc.Connectivity = connectivity;
cc.ImageSize = sz;
cc.NumObjects = 1;
cc.PixelIdxList = {ind};
