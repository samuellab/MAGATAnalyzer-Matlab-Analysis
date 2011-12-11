function [kernel, amp, sigma, offset] = findGaussianKernel(xdata, convdata, x)
%function [kernel, amp, sigma, offset] = findGaussianKernel(xdata, convdata)
%
%finds the gaussian kernel that minimizes the square error between
%conv(xdata, kernel) and convdata
%amp = area of gaussian
%sigma = standard deviation of gaussian
%offset = location of maximum point relative to center of kernel

existsAndDefault('x',[mean(convdata)/mean(xdata), 1, 0]);
x = lsqcurvefit(@convWithKernel,x,xdata,convdata);
kernel = gaussWithOffset(x(1), x(2), x(3));
amp = x(1);
sigma = x(2);
offset = x(3);

function yout = convWithKernel(x, xdata)
k = gaussWithOffset(x(1), x(2), x(3));
yout = conv(xdata, k, 'same');
