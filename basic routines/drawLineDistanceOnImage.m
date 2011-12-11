function [im,valueim] = drawLineDistanceOnImage(im, pt1, pt2, valueim, value)
%sets all points that are on a line running
%perpendicular to the line pt1-pt2 that passes through the segment 
%between pt1&pt2 to the value of the distance to that line segment
%
%note that the x,y point of im is found in im(y,x)
%while pt1,pt2 should have form [x y]
%
%pt1,pt2 can be a vector of points, in which case we will draw the minimum
%of all distances
%
%optional parameter valueim, value
%if valueim is non-empty, then we put value(j) in valueim wherever the
%distance between the point and the segment btw pt1(j) & pt2(j) is the
%minimum distance between the point and all other segments
%
%to speed things up, we only consider points within a reasonable range of
%the line segment

if (size(pt1,1) ~= 2)
    pt1 = pt1';
end
if (size(pt2,1) ~= 2)
    pt2 = pt2';
end

v = pt2-pt1;
l = sqrt(sum(v.^2));
v = single(v./[l;l]);
n = single([-v(2,:);v(1,:)]);

if (existsAndDefault('valueim', []))
    existsAndDefault('value', 1:length(pt1));
    if (length(value) == 1)
        value = repmat(value, 1, length(pt1));
    end
end


[x,y] = meshgrid(1:size(im,2),1:size(im,1));

d = (sqrt(sum((pt2-pt1).^2, 1))); %distance between each set of points

xlim = round(sort([pt1(1,:);pt2(1,:)],1) + 2*[-abs(d.*n(1,:));abs(d.*n(1,:))]);
ylim = round(sort([pt1(2,:);pt2(2,:)],1) + 2*[-abs(d.*n(2,:));abs(d.*n(2,:))]);

xlim(xlim<1) = 1;
xlim(xlim>size(x,2)) = size(x,2);
ylim(ylim<1) = 1;
ylim(ylim>size(x,1)) = size(x,1);



for j = 1:size(pt1,2)

     %try
        xx = x(ylim(1,j):ylim(2,j), xlim(1,j):xlim(2,j));
        yy = y(ylim(1,j):ylim(2,j), xlim(1,j):xlim(2,j));
     %catch me
      %   disp(me.getReport);
      %   size(x)
       %  xlim(1,j)
       %  xlim(2,j)
        % ylim(1,j)
        
         %ylim(2,j)
     %end
    sz = size(xx);
    
    xx = reshape(xx,1,[]);
    yy = reshape(yy,1,[]);
    

    frompt1 = single([xx-pt1(1,j);yy-pt1(2,j)]);
    frompt2 = single([xx-pt2(1,j);yy-pt2(2,j)]);
    
    vv = repmat(v(:,j), 1, length(xx));
    nn = repmat(n(:,j), 1, length(xx));
    dist = abs(dot(frompt1,nn));
    
    online = dot(vv, frompt1) >= 0 & dot(vv,frompt2) <= 0;    
    [I,J] = ind2sub(sz, find(online));
    inds = sub2ind(size(im), I+ylim(1,j)-1, J + xlim(1,j) - 1);
    
    if (~isempty(valueim))
        valueim(inds(dist(online) < im(inds))) = value(j);
    end
    im(inds) = min(im(inds), dist(online));
    
end

