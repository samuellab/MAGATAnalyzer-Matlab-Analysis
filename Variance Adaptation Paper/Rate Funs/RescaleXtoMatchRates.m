clear rates_err

nBins = 5;
varL = or42a_R2.var;
varH = or42a_R8.var;

lxs = @(alpha, varstruct) linspace(-3*std(alpha*varstruct.noturn.x_conv), 3*std(alpha*varstruct.noturn.x_conv), nBins);

alphas = linspace(1, 2, 30);

lxL = lxs(1, varL);
lxH = lxs(1, varH);
[trL, trL_eb] = turnRate(1, varL, nBins);
[trH, trH_eb] = turnRate(1, varH, nBins);

for i=1:length(alphas)
    tic
    a = alphas(i);
    
    lxL_s{i} = lxs(a, varL);
    trL_s{i} = turnRate(a, varL, nBins);
    
    if(lxL_s{i}(end)>lxH(end))
        rates_err(i) = immse(trH, interp1(lxL_s{i}, trL_s{i}, lxH));
    else
        rates_err(i) = immse(trL_s{i}, interp1(lxH, trH, lxL_s{i}));
    end
   toc 
end


% op = optimoptions('fmincon');
% op.Display = 'off';
% op.Algorithm = 'quasi-newton';
% problem.options = op;
% problem.solver = 'fmincon';
% problem.objective = dist;
% problem.x0 = 1;
% problem.lb = .6;
% problem.ub = 1.4;


function [rate, rate_eb] = turnRate(alpha, varstruct, nBins)

deltaT = median(diff(varstruct.fullensemble.eti));

xt = alpha * varstruct.turn.x_conv;
xr = alpha * varstruct.noturn.x_conv;

lx = linspace(-3*std(xr), 3*std(xr), nBins);

nt = histcounts(xt, binEdgesFromCenters(lx));
nr = histcounts(xr, binEdgesFromCenters(lx));

rate = 1/deltaT * nt ./ nr;
rate_eb = 1/deltaT * sqrt(nt) ./ nr;

op = optimoptions('fminunc');
op.Display = 'off';
op.Algorithm = 'quasi-newton';

nlogp = @(p) -sum(polyval(p, xt)) + deltaT*sum(exp(polyval(p,xr)));
p0 = [mean(xt)/var(xr) log(length(xt)/(length(xr)*deltaT))];

fP = fminunc(nlogp, p0, op);
rate = exp(polyval(fP, lx));

end

