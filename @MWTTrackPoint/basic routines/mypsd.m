function [ps, k,fullps] = mypsd (data, npoints)
%function [ps, k, fullps] = mypsd (data, npoints)
%

k = (0:(npoints-1))/npoints;


paddata = zeros([1 ceil(length(data)/npoints)*npoints]);
paddata(1:length(data)) = data;

paddata = reshape (paddata, npoints, []);

ffdata = fft(paddata,[],1) / (npoints/2);

ps = mean(abs(ffdata).^2, 2);
fullps = ps;
k = k(1:(npoints/2));
ps = ps(1:(npoints/2));