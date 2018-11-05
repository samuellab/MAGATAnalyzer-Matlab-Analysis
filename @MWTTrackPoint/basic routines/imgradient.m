function [xderiv,yderiv, tderiv] = imgradient (im, sigma)
%function [xderiv,yderiv] = imgradient (im, sigma)
if (size(im, 3) > 1)
    xderiv = zeros(size(im));
    yderiv = xderiv;
    for j = 1:size(im,3)
        [xderiv(:,:,j), yderiv(:,:,j)] = imgradient(im(:,:,j), sigma);
    end
    if (nargout > 2)
        tderiv = zeros(size(xderiv));
        for j = 1:size(im,1)
            tderiv(j,:,:) = imgradient(squeeze(im(j,:,:)), sigma);
        end
    end
    return;
end
gk = gaussKernel(sigma);
dg = dgausskernel(sigma);
padsize = floor(length(gk)/2);
padim = padarray(im, [padsize padsize], 'replicate');

xderiv = conv2(gk, dg, double(padim), 'valid');
if (nargout > 1)
    yderiv = conv2(dg, gk, double(padim), 'valid');
end