function [x, c, ctr] = fitFourierToImage (im, n, x0)
%c = fitLegendrePolynomialToImage (im, n)
im = double(im);
t = linspace(0,2*pi,100);
p = zeros(2*n+1,length(t));
p(1,:) = 1;
for j = 1:n
%    P = legendre(j-1,cos(t/2));  
    p(2*j,:) = cos(j*t);
    p(2*j+1,:) = sin(j*t);
end
op = optimset('fminunc');
op.DiffMinChange = 0.25;
op.Display = 'off';
op.LargeScale = 'off';
if (nargin < 3)
    c = zeros(1, 2*n+1);
    xx = repmat(1:size(im,2), size(im,1), 1);
    yy = repmat((1:size(im,1))', 1, size(im,2));
    x = [sum(xx(:).*im(:)) sum(yy(:).*im(:))]/sum(im(:));
    c(1) = sqrt((sum((xx(:)-x(1)).^2.*im(:))+ sum((yy(:)-x(2)).^2.*im(:)))/sum(im(:)));
    c(2:end) = c(1)*0.1*randn(size(c(2:end)));%c(1)*0.1;
    x0 = [x c];
    
    for s = [min(size(im))/6, min(size(im))/12, 1, 0.1]
        myfun = @(q) energyfun(q, blurim(im,s), p, t);
        x0 = fminunc(myfun, x0,op);
    end
end
scale = 3;
im2 = imresize(im,scale);
myfun = @(q) energyfun(q, im2, p, t);
x1 = fminunc(myfun, x0*scale,op);
x = x1(1:2)/scale;
c = x1(3:end)/scale;
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
% imagesc(im); axis equal; colorbar vert; hold on; plot(ctrx, ctry, 'r-'); hold off; pause(0.01);
en = sum((im(bwim)-mean(im(bwim))).^2) + 2*sum((im(~bwim) - mean(im(~bwim))).^2); 

end