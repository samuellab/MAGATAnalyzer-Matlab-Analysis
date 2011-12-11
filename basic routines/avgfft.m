function [f,v] = avgfft (data, blocksize)
%function [f,v] = avgfft (data, blocksize)


f = zeros([1 blocksize]);
%n = 0;
for j = 1:blocksize:length(data)
    nextind = min(j + blocksize - 1, length(data));
    f = f + fft(data(j:nextind), blocksize);
 %   n = n+1;
end

f = f / length(data);
v = ifftshift(linspace(-1,1,blocksize)/2);
 