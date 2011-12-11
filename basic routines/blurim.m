function im = blurim (im, sigma)
%function im = blurim (im, sigma)
gk = gaussKernel(sigma);
padsize = floor(length(gk)/2);
padim = padarray(im, [padsize padsize], 'replicate');

im = conv2(gk, gk, padim, 'valid');