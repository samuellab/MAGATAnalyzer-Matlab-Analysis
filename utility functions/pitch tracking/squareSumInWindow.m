function ssf = squareSumInWindow (x,minshift, maxshift, winsize, gaussWindow) 
%function ssf = squareSumInWindow (x,minshift, maxshift, winsize, gaussWindow) 
%
%written for 1D data right now

existsAndDefault('gaussWindow', false);
minshift = floor(minshift);
maxshift = ceil(maxshift);
nrows = maxshift - minshift + 1;
n = length(x);
ssf = zeros(nrows, n);
for j = 1:nrows
    inds = (1:(n-j-minshift + 1));
    ssf(j,inds) = (x(inds)).^2+(x(inds + minshift + j -1)).^2;
end

if gaussWindow
    k = gaussKernel(winsize);
else
    k = ones(1,ceil(winsize))/ceil(winsize);
end

for j = 1:size(ssf, 1)
    ssf(j,:) = conv2(ssf(j,:), k, 'same');
end


 