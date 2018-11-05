function dim = distanceFromOrientedContourImage (cpts, imOrigin, sizeofim)
%function dim = distanceFromOrientedContourImage (cpts, imOrigin, sizeofim)
%
%distance from edge of contour to image points; dim > 0 if points are
%outside; dim < 0 for points inside

if (~all(cpts(:,1) == cpts(:,2)))
    cpts = cpts(:, [1:end 1]);
end
cpts(1,:) = cpts(1,:)-imOrigin(1) + 1;
cpts(2,:) = cpts(2,:)-imOrigin(2) + 1;

%pm = poly2mask(cpts(1,:), cpts(2,:), sizeofim(1), sizeofim(2));

[x,y] = meshgrid((1:sizeofim(2)), (1:sizeofim(1)));
allpts = [x(:) y(:)]';
dim = abs(distanceFromPointToLineSegment(cpts(:,1), cpts(:,2), allpts));

for j = 2:(size(cpts,2)-1)
    d = abs(distanceFromPointToLineSegment(cpts(:,j), cpts(:,j+1), allpts));
    dim(d < dim) = d(d < dim);
end
dim = reshape(dim, size(pm));
dim(pm) = -dim(pm);
