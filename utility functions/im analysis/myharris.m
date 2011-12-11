function [cim, cpts] = myharris (im, sig)
%function [cim, cpts] = myharris (im, sig)
%
%if sigma is a vector, finds corners iteratively, using largest sigma first

sig = sort(sig, 'descend');

cimbwold = true(size(im));

for j = 1:length(sig)
    sigma = sig(j);

    [x,y] = imgradient(im, sigma);
    ix2 = blurim(x.^2, sigma);
    iy2 = blurim(y.^2, sigma);
    ixy = blurim(x.*y, sigma);

    cim = (ix2.*iy2 - ixy.^2)./(ix2 + iy2 + eps);

    if nargout <= 1
        return;
    end

    [idx,c] = kmeans(double(cim(:)), 2, 'start', double([min(cim(:)); max(cim(:))]));
    if (c(2) > c(1))
        high = 2;
    else
        high = 1;
    end
    cimbw = cimbwold & reshape(idx == high, size(cim));
    cimbwold = cimbw;
end
stats = regionprops(cimbw, abs(cim), 'WeightedCentroid');
pts = [stats.WeightedCentroid];
cpts = [pts(1:2:end); pts(2:2:end)];    