function c = shimazakiCostFunction(turnTimes, w)
%function c = shimazakiCostFunction(sortedTurnTimes, w)
%c = N/w + 2/w sum(i < j) exp(-(ti - tj)^2/
%
tt = repmat(turnTimes(:), [1 length(turnTimes)]);

dtsw = (tt - tt').^2/(4*w^2);
dtsw = dtsw(:);
c = 1/w*(length(turnTimes)*2*sqrt(2) + sum(exp(-dtsw) - 2*sqrt(2)*exp(-2*dtsw)));
