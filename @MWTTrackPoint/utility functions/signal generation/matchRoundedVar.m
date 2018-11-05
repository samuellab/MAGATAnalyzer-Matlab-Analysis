function normSigma = matchRoundedVar(targetSigma)
%function normSigma = matchRoundedVar(targetSigma)
%gives a new value of sigma s.t. round(randn([1 N]) * normSigma) has the
%standard deviation targetSigma
%
%note that despite function name, you give target std dev, not variance

if (length(targetSigma) > 1)
    normSigma = targetSigma;
    for j = 1:length(targetSigma)
        normSigma(j) = matchRoundedVar(targetSigma(j));
    end
    return;
end

if (2*normcdf(-1.5, 0, targetSigma) < 1e-4) %probability of a number other than 0,+/-1 < 1e-4, so variance is = pval of 0.5 | sigma
    normSigma = -0.5./norminv(.5 * targetSigma.^2, 0, 1);
    return;
end

%maxreps = 100;
nsamples = 1e7;
%normSigma = targetSigma;
rd = randn([1 nsamples]);
%convergence = 1e-3;

normSigma = fminbnd(@(s) (std(round(rd*s))-targetSigma).^2, targetSigma/10, targetSigma * 2);

% for j = 1:maxreps
%     s = std(round(rd*normSigma));
%     if (s == 0)
%         normSigma = normSigma * 2;
%     else
%         normSigma = normSigma * targetSigma/s;
%     end
%     if (abs(s/targetSigma - 1) < convergence)
%         break;
%     end
% end