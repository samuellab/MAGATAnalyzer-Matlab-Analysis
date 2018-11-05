function pval = binomialpval (m,n,p)
%function pval = binomialpval (m,n,p)
%
%if the null hypothesis is that the p fraction of the time, the result is
%true, then
%pval is the probability at least m out of n trials are true

sigma = sqrt(n*p.*(1-p));
mn = n*p;
pval = 1/2*(erfc((m-mn)./(sqrt(2)*sigma)));
