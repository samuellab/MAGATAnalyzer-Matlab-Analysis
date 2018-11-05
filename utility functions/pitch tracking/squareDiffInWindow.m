function sdf = squareDiffInWindow (x,minshift, maxshift, winsize, gaussWindow) 
%function sdf = squareDiffInWindow (x,minshift, maxshift, winsize, gaussWindow) 
%
%written for 1D data right now

existsAndDefault('gaussWindow', false);
minshift = floor(minshift);
maxshift = ceil(maxshift);
nrows = maxshift - minshift + 1;
n = length(x);
sdf = zeros(nrows, n);
for j = 1:nrows
    inds = (1:(n-j-minshift + 1));
    sdf(j,inds) = (x(inds)-x(inds + minshift + j -1)).^2;
end

if gaussWindow
    k = gaussKernel(winsize);
else
    k = ones(1,ceil(winsize))/ceil(winsize);
end

for j = 1:size(sdf, 1)
    sdf(j,:) = conv2(sdf(j,:), k, 'same');
end


 