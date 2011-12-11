function [cunb,c, np] = xcorrVec (a, b)
%a,b are m x N matrices (or N x m), where m is the dimension
%and N is the number of vectors in the time series and m < N
%for instance, a sequence of 1000 x,y points would be a 2x1000 vector
%
%%if a and b are different lengths, the shorter vector is 
%zero padded (at the end?)
%
%we define the cross-correlation to be
%xc(T) = <dot(a(t),b(t-T))>; -N<T<N and <> denotes the average
%c is the return 1x(2N-1) vector c(j) = xc(j-N);
%
%the autocorrelation for T = 0 to N is given by c = xcorrVec(a,a); ac = c(N:end);
%cunb is the unbiased (normalized) correlation
%c is the unnormalized correlation
%np is the number of points used to calculate c or cunb;  cunb = c./np
na = length(a);
nb = length(b);
n = max(na, nb);
c = zeros([1 2*n-1]);

if (size(a,1) > size(a,2))
    a = a';
end
if (size(b,1) > size(b,2))
    b = b';
end

for k = 1:size(a,1)
    c = c+xcorr(a(k,:),b(k,:));
end
np = ([(1:n) ((n-1):-1:1)]);
cunb = c./np;

