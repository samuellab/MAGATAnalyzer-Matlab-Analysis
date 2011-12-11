function im = drawLineOnImage(im, pt1, pt2, thickness, value)
%sets all points within (thickness/2) pixels on the line between pt1 and pt2  
%to value
%
%note that the x,y point of im is found in im(y,x)
%while pt1,pt2 should have form [x y]
%
%pt1,pt2 can be a vector of points, in which case we will draw lines
%between all sets of pt1 and pt2
%if value is also a vector, then we draw value(k) between pt1(k) & pt2(k)


if (size(pt1,1) ~= 2)
    pt1 = pt1';
end
if (size(pt2,1) ~= 2)
    pt2 = pt2';
end

if (length(value) == 1)
    value = repmat(value, 1, length(pt1));
end

v = pt2-pt1;
l = sqrt(sum(v.^2));
v = single(v./[l;l]);
n = single([-v(2,:);v(1,:)]);


[x,y] = meshgrid(1:size(im,2),1:size(im,1));

xlim = round(sort([pt1(1,:);pt2(1,:)],1) + repmat([-thickness;thickness],1,length(pt1)));
ylim = round(sort([pt1(2,:);pt2(2,:)],1)+ repmat([-thickness;thickness],1,length(pt1)));

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
    
    online = dot(vv, frompt1) > 0 & dot(vv,frompt2) < 0 & abs(dot(frompt1,nn)) < thickness/2;
    
    [I,J] = ind2sub(sz, find(online));
    inds = sub2ind(size(im), I+ylim(1,j)-1, J + xlim(1,j) - 1);
    
    
    im(inds) = value(j);

end

