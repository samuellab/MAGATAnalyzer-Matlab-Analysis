function ms = skewNormalMeanSq(u, s, a)
%
%mean and standard deviation of skew-normal distribution with location u, scale s, and skew-parameter a
%http://en.wikipedia.org/wiki/Skew_normal_distribution
[mu, sigma] = skewNormalMean(u, s, a);
ms = mu.^2 + sigma.^2;
end

