function [freq, phase, ampl, yfit] = fitWave(xd, ydata, ramptype, varargin)
% function [freq, phase, ampl] = fitWave(xdata, ydata, ramptype, varargin)
%
% finds parameters of a wave that best describes the data
% ramptype is one of 'square', 'triangle', 'sine'
%
% varargin: 'fixedPeriod', t - if the period is known, input it here.

fixedPeriod = [];
varargin = assignApplicable(varargin);

[xd, I] = sort(xd);
ydata = ydata(I);

dx = diff([(2*xd(1) - xd(2)) xd]);
ydm = sum(ydata.*dx)/sum(dx);

yd = ydata - ydm;

switch (lower(ramptype))
    case 'square'
        yfn = @(x, xdata) sign(sin(x(1)*xdata + x(2)));
    case 'sine'
        yfn = @(x, xdata) sin(x(1)*xdata + x(2));
    case 'triangle'
%         if (exist('sawtooth', 'file'))
%             yfn = @(x, xdata)sawtooth(x(1)*xdata + x(2), 0.5);
%         else
%             warning ('signal processing toolbox not installed? using sine for triangle wave');
            yfn = @(x, xdata) mytriangle(x(1)*xdata + x(2));
        %end
    otherwise
        yfn = @(x, xdata) sin(x(1)*xdata + x(2));
end

%testfun = @(x) -sum(yfn(x,xd).*yd.*dx);

xi = linspace(xd(1), xd(end), length(xd));
yi = interp1(xd, yd, xi);
fs = 2*pi/mean(diff(xi));
%adapted from matlab help doc

NFFT = 2^nextpow2(length(yi)); % Next power of 2 from length of y
Y = fft(yi,NFFT)/length(yi);
f = fs/2*linspace(0,1,NFFT/2+1);
[~,I] = max(abs(Y(1:(NFFT/2+1))));
x(1) = f(I); %initial guess for the frequency - highest power in spectrum
%initial guess for phase - find minimum to within 1%
pp = linspace(0,2*pi,100); tf = zeros(size(pp));
testfun = @(x) -sum(yfn(x,xd).*yd.*dx);
if (~isempty(fixedPeriod))
    x(1) = 2*pi/fixedPeriod;
end
% for j = 1:length(pp), tf(j) = testfun([x(1), pp(j)]); end
% [~,II] = min(tf);
% x(2) = pp(II)
% 
x(2) = fminbnd(@(y) testfun([x(1) y]), 0, 2*pi);
xinit = x;
if (isempty(fixedPeriod))
    x = fminsearch(testfun, xinit); %keeping xinit around for debugging purposes - note fminsearch succeeds where fmincon fminunc fail
end
freq = x(1);
phase = x(2);
ampl = regress(yd',yfn(x, xd)');
yfit = ampl * yfn(x,xd) + ydm;



            