
function [btdstruct, RatesinCycle] = findRates(btdstruct, opstruct, T, tau)
% function [btdstruct, RatesinCycle] = findRates(btdstruct, opstruct, T, tau)

% Finds high and low rate functions for time bins of duration T within each cycle (defined by btdstruct.var(1).period),
% Excludes the first tau seconds. Does it for first and second 10min separately

if (~isfield(btdstruct, 'var'))
    btdstruct = BehaviorTriggeredData.prepVarianceSwitchingAnalysis_Gepner(btdstruct, opstruct);
end
var = btdstruct.var;

if(~isfield(btdstruct, 'varops'))
    fn = {'stim', 'numLxBins', 'adaptationTime', 'Q_alpha', 'numTimeBins', 'timeRange'};
    for j = 1:length(fn)
        btdstruct.varops.(fn{j}) = opstruct.(fn{j});
    end
end

nstim = size(var, 2);
period = var(1).period;

for j=1:nstim
    deltaT = median(diff(var(j).fullensemble.eti));
    tton = var(j).turn.ton;
    teti = var(j).turn.eti;
    rton = var(j).noturn.ton;
    reti = var(j).noturn.eti;
    
    
    n = period - tau;
    n = round(n/T);
    
    for i=1:n
        t1 = tau + (i-1)*T;
        t2 = tau + i*T;
        
        rH = rton>=t1 & rton<=t2 & reti >= min(opstruct.timeRange) & reti < max(opstruct.timeRange);
        tH = tton>=t1 & tton<=t2 & teti >= min(opstruct.timeRange) & teti < max(opstruct.timeRange);
        xt = var(j).turn.x_conv(tH);
        xr = var(j).noturn.x_conv(rH);
        lx = linspace(-3*std(xr), 3*std(xr), opstruct.numLxBins);
        RatesinCycle{i}(j) = turnRates (xt, xr, deltaT, lx);
        
        
    end
    
    
end

btdstruct.RatesinCycle = RatesinCycle;

end

function rf = turnRates (xt, xr, deltaT, lx)

op = optimoptions('fminunc');
op.Display = 'off';
op.Algorithm = 'quasi-newton';

existsAndDefault('lx', linspace(-3*std(xr), 3*std(xr), 20));

nt = histcounts(xt, binEdgesFromCenters(lx));
nr = histcounts(xr, binEdgesFromCenters(lx));

rf.lx = lx;
rf.rate = 1/deltaT * nt ./ nr;
rf.rate_eb = 1/deltaT * sqrt(nt) ./ nr;

nlogp = @(p) -sum(polyval(p, xt)) + deltaT*sum(exp(polyval(p,xr)));
p0 = [mean(xt)/var(xr) log(length(xt)/(length(xr)*deltaT))];
rf.rateFitLin = fminunc(nlogp, p0, op);
p0 = [0 rf.rateFitLin];
rf.rateFitQuad = fminunc(nlogp, p0, op);

end 





