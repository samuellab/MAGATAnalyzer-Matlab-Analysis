function pval = binomialpval_exact (m,n,p)
%function pval = binomialpval_exact (m,n,p)
%
%if the null hypothesis is that the p fraction of the time, the result is
%true, then
%pval is the probability at least m out of n trials are true


pval = zeros(size(m));
if (length(p) == 1)
    p = repmat(p,size(m));
end
for j = 1:length(m)
    for k = m(j):n(j)
        pval(j) = pval(j) + nchoosek(n(j),k)*p(j)^(k)*(1-p(j))^(n(j)-k);
    end
end