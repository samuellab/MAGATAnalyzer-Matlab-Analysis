function im = traceBoundarySubPixel (cpts, imOrigin, sizeofim)
%function im = traceBoundarySubPixel (cpts, imOrigin, sizeofim)
%
%finds all pixels that intersect the boundary line, then labels each of
%these pixels with the fraction (between 0.0 and 1.0) of the pixel
%contained within the contour.  
%
%im(i,j) represents the pixel spanning the rectangle (i + imOrigin(1),j+imOrigin(2)) and (i+imOrigin(1) + 1,j+imOrigin(2)+1)

if (~all(cpts(:,1) == cpts(:,2)))
    cpts = cpts(:, [1:end 1]);
end
x = cpts(1,:)-imOrigin(1) + 1;
y = cpts(2,:)-imOrigin(2) + 1;

%close the contour
if (x(end) ~= x(1) || y(end) ~= y(1))
    x(end+1) = x(1);
    y(end+1) = y(1);
end

im = double(poly2mask(x,y,sizeofim(1), sizeofim(2)));
%fails at corners!
im = markBoundary (im, x,y);

%do corners separately
xx = x(2:end);
yy = y(2:end);
for j = 1:length(xx)
    if (xx(j) < 1 || floor(xx(j)) > sizeofim(2) || yy(j) < 1 || yy(j) > sizeofim(1))
        continue;
    end
    ind1 = mod(j-2, length(xx))+1;
    ind2 = mod(j,length(xx)) + 1;
    im(floor(yy(j)),floor(xx(j))) = divideSquareByCorner(floor([xx(j);yy(j)]), floor([xx(j);yy(j)])+[1;1], [xx(ind1);yy(ind1)], [xx(j);yy(j)], [xx(ind2);yy(ind2)]);
end
if any(size(im) ~= sizeofim)
    error ('something bad happened');
end

%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here


end

