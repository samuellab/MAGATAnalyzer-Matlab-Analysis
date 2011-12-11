function [c,com] = imCov(im, com)
%function [c,com] = imCov(im, com)
%
%c is the covariance matrix defined by 
%c(i,j) mean(Im(k) x(i,k) * x(j,k))
%where x(i,k) = pos(i,k) - com(i)
%
%if com is not passed in, it is calculated by
%com = mean(Im(k) pos(:,k))
im = double(im);
[x,y] = meshgrid(1:size(im,2), 1:size(im,1));
if (~exist('com', 'var'))
    com(1) = sum(x(:) .* im(:))/sum(im(:));
    com(2) = sum(y(:) .* im(:))/sum(im(:));
end;

x = x(:) - com(1);
y = y(:) - com(2);

c = zeros(2,2);
c(1,1) = sum(im(:) .* x.^2);
c(1,2) = sum(im(:) .* x .* y);
c(2,1) = c(1,2);
c(2,2) = sum(im(:) .* y.^2);
c = c/sum(im(:));
