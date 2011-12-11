function [xderiv,yderiv] = imgradient (im, sigma)
%function [xderiv,yderiv] = imgradient (im, sigma)
gk = gaussKernel(sigma);
dg = dgausskernel(sigma);
padsize = floor(length(gk)/2);
padim = padarray(im, [padsize padsize], 'replicate');

xderiv = conv2(gk, dg, double(padim), 'valid');
yderiv = conv2(dg, gk, double(padim), 'valid');