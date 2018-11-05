function dkl = DKL_gauss_est(x0, x1 )
%function dkl = DKL_gauss_est(x1, x2 )
%
%estimates the kl divergence between set x1 and set x2 using the mean and
%covariances only (e.g. approximates the distributions as gaussian)
% dkl = dkl(x1||x2)
%
%size(x1) and size(x2) should be MxD and NxD where M,N are the number of
%samples and D is the dimension

if (size(x0,2) ~= size(x1,2))
    error ('second dimension must be the same for x1,x2');
end
v0 = cov(x0,1);
u0 = mean(x0);
v1 = cov(x1,1);
u1 = mean(x1);


dkl = 0.5*(trace(v1\v0) + (u1-u0)*(v1\(u1-u0)') - size(x0,2) + log(det(v1))-log(det(v0)));

end

