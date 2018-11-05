function [pval, pcommon] = binomialpval_twogroups (na, ka, nb, kb)
%function [pval, pcommon] = binomialpval_twogroups (na, ka, nb, kb)
%
%group a has na trials with ka successes
%group b has nb trials with kb successes
%assume (ka/na) > (kb/nb)
%
%pval is the probability that na trials with a pcommon probability of
%success produces AT LEAST ka successes AND that nb trials with pcommon
%probability of success produces AT MOST kb successes.  
if (ka > na || kb > nb)
    error ('BPV:badval', 'number of successes cannot be greater than number of trials');
end

if (ka/na < kb/nb)
    ktemp = ka;
    ntemp = na;
    ka = kb;
    na = nb;
    kb = ktemp;
    nb = ntemp;
end

pvfun = @(p) -((1-binocdf(ka-1, na, p))*binocdf(kb,nb,p));

%pcommon = (ka+kb)/(na+nb);

[pcommon,pval] = fminbnd(pvfun, 0, 1);
pval = -pval;