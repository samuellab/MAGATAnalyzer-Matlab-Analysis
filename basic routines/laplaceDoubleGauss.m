function k = laplaceDoubleGauss (sigma1, sigma2)
%function k = laplaceDoubleGauss (sigma1, sigma2)
%creates a 2D laplacian kernel by subtracting a regional gaussian (sigma2)
%from a local gaussian (sigma1).  

sigma = max(sigma1,sigma2);

x = floor(-3*sigma):ceil(3*sigma);
g = exp(-x.^2/(2 * sigma1));
g = g'*g;
g1 = g / (sum(sum(g)));

g = exp(-x.^2/(2 * sigma2));
g = g'*g;
g2 = g / (sum(sum(g)));

k = g1 - g2;
