function [mu, sigma] = skewNormalMean(u, s, a)
%function [mu, sigma] = skewNormalMean(u, s, a)
%mean and standard deviation of skew-normal distribution with location u, scale s, and skew-parameter a
%http://en.wikipedia.org/wiki/Skew_normal_distribution
delta = a./sqrt(1+a.^2);
mu = u + delta.*s.*sqrt(2/pi);
sigma = s.*sqrt(1-2.*delta.^2/pi);
end

