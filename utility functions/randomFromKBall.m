function x = randomFromKBall (k, npts)
%function x = randomFromKBall (k, npts)
%
%picks random points uniformly from the ball in k dimensions

x = randn(k, npts);
r = sqrt(sum(x.^2, 1));
pf = rand(1, npts).^(1/k)./r;
for j = 1:k
    x(j,:) = x(j,:).*pf;
end
