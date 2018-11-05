function acf = autocorrInWindow (x,minshift, maxshift, winsize, gaussWindow) 
%function acf = autocorrInWindow (x,minshift, maxshift, winsize, gaussWindow) 
%
%written for 1D data right now

existsAndDefault('gaussWindow', false);
minshift = floor(minshift);
maxshift = ceil(maxshift);
nrows = maxshift - minshift + 1;
n = length(x);
acf = zeros(nrows, n);
for j = 1:nrows
    inds = (1:(n-j-minshift + 1));
    acf(j,inds) = x(inds).*x(inds + minshift + j -1);
end

if gaussWindow
    k = gaussKernel(winsize);
else
    k = ones(1,ceil(winsize))/ceil(winsize);
end

for j = 1:size(acf, 1)
    acf(j,:) = conv2(acf(j,:), k, 'same');
end


 