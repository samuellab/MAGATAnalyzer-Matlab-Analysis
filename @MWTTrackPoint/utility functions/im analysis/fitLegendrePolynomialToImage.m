function [x, c, ctr] = fitLegendrePolynomialToImage (im, n)
%c = fitLegendrePolynomialToImage (im, n)
im = double(im);
t = linspace(0,2*pi,100);
p = zeros(n,length(t));
for j = 1:n
    P = legendre(j-1,cos(t/2));
    p(j,:) = P(1,:);
end
c = zeros(1, n);
xx = repmat(1:size(im,2), size(im,1), 1);
yy = repmat((1:size(im,1))', 1, size(im,2));
x = [sum(xx(:).*im(:)) sum(yy(:).*im(:))]/sum(im(:));
c(1) = 2*sqrt((sum((xx(:)-x(1)).^2.*im(:))+ sum((yy(:)-x(2)).^2.*im(:)))/sum(im(:)));
c(2:end) = c(1)*0.1;
x0 = [x c];
myfun = @(q) energyfun(q, im, p, t);
myfun(x0)
x1 = fminsearch(myfun, x0);
x = x1(1:2);
c = x1(3:end);
[ctrx,ctry] = lpToCtr(x,c,p,t);
ctr = [ctrx;ctry];
    
end

function en = energyfun (q, im, p, t)
    [ctrx,ctry] = lpToCtr(q(1:2), q(3:end), p, t);
    en = energyInContour(im, ctrx, ctry);
end    
    
function [ctrx,ctry] = lpToCtr(x,c,p,t)
    r = c*p;
    ctrx = x(1) + r.*cos(t);
    ctry = x(2) + r.*sin(t);
end

function en = energyInContour (im, ctrx,ctry)
 %   bwim = traceBoundarySubPixel([ctrx;ctry],[1 1], size(im));

  %  en = sum((im(:) - mean(im(:).*bwim(:))).^2.*bwim(:)) ...
  %      + sum((im(:) - mean(im(:).*(1-bwim(:)))).^2.*(1-bwim(:)));
%    imagesc(bwim); axis equal; colorbar vert; pause(0.01);

bwim = poly2mask(ctrx, ctry, size(im,1), size(im,2));
 imagesc(bwim); axis equal; colorbar vert; pause(0.01);
en = sum((im(bwim)-mean(im(bwim))).^2) + sum((im(~bwim) - mean(im(~bwim))).^2); 

end