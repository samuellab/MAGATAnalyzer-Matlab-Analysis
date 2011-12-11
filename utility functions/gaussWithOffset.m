function kernel = gaussWithOffset(amp, sigma, offset)
%function kernel = gaussWithOffset(amp, sigma, offset)

SIGMA = 20;
k = gaussKernel(SIGMA);
kernel = zeros([1, ceil(6*sigma + 1 + abs(2*offset))]);
kinds = linspace(-(length(k)-1)/2, (length(k)-1)/2, length(k))*sigma/SIGMA;
kerninds = linspace(-(length(kernel)-1)/2, (length(kernel)-1)/2, length(kernel)) - offset;
%{
min(kinds)
max(kinds)
min(kerninds)
max(kerninds)
%}
kernel = interp1(kinds, k, kerninds, 'linear', 0);
kernel = amp*kernel/sum(kernel);