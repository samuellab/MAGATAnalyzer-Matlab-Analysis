function problemDescription = simulateVarSwitch (btd, opstruct, oldpd, timeField, trange, staticParams, nlarvae)
% problemDescription = simulateVarSwitch (btd, opstruct, timeField, polynomialDegree, trange, staticParams, temporalParams, nlarvae)
%

if (nargin < 1)
    opstruct.stim(1).gqname = 'led1ValDiff';
    opstruct.stim(1).iswn = true;
    opstruct.stim(1).ramptype = 'square';
    opstruct.stim(1).switch = true;
    opstruct.stim(1).period = [];
    opstruct.taxis = 60:.05:1200;
    opstruct.alphaLow = 0.7;
    opstruct.alphaHigh = 1.3;
    problemDescription = opstruct;
    disp ('ramp types are square, triangle, sine, and constant');
    return;
end
fn = {'ratefun', 'gradlogratefun', 'hesslogratefun', 'params_0', 'temporalratemod'};
for j = 1:length(fn)
    pd.(fn{j}) = oldpd.(fn{j});
end
existsAndDefault('trange', oldpd.timeRange);

for j = find(~[opstruct.stim.iswn] & [opstruct.stim.switch])
    disp (['stim ' - num2str(j) ' is not switching or noisy - so I am ignoring it']);
end

varswitch = [opstruct.stim.iswn] & [opstruct.stim.switch];
switchstim = 0;
switch (nnz(varswitch))
    case 0
        if (any([opstruct.stim.switch]))
            switchstim = find([opstruct.stim.switch], 1, 'first');
            gq_switch = opstruct.stim(switchstim).gqname;
        end
    case 1
        switchstim = find(varswitch);
        gq_switch = [opstruct.stim(switchstim).gqname '_var'];
    otherwise
        warning ('I don''t know how to handle multiple variances switching simultaneously');
        return;
end
if (strcmpi(opstruct.stim(switchstim).ramptype, 'square'))
    gqlow = [gq_switch, '_low'];
else
    gqlow = [sq_switch, '_rising'];
end

ndim = nnz([opstruct.stim.iswn]);

if(switchstim==0)
    xnames = cellfun(@(s) [s 'TurnLin'], {opstruct.stim([opstruct.stim.iswn]).gqname}, 'UniformOutput', false);
else
    xnames = cellfun(@(s) [s 'TurnLin'], {opstruct.stim([opstruct.stim.switch]).gqname}, 'UniformOutput', false);
end

xnames = [xnames gqlow];

dt = min(0.05, median(diff(opstruct.taxis)));

pd.rateTx = min(opstruct.taxis):dt:max(opstruct.taxis);

if(isfield(oldpd, 'alpha_eti'))
    
    tx_eti = oldpd.tx_eti;
    alpha_eti = oldpd.alpha_eti;
    alpha_eti = interp1(tx_eti, alpha_eti, pd.rateTx);
    %     ratefun = @(xdata, tdata) pd.temporalratemod(temporalParams, tdata) .* pd.ratefun(staticParams, alpha_eti.*xdata{1});
    ratefun = @(xdata, tdata) pd.ratefun(staticParams, alpha_eti.*xdata{1});
else
    
    alphafun = @(xd) opstruct.alphaHigh*double(logical(xd)) + opstruct.alphaLow*double(~logical(xd));
    %     ratefun = @(xdata, tdata) pd.temporalratemod(temporalParams, tdata) .* pd.ratefun(staticParams, alphafun(xdata{end}).*xdata{1});
    ratefun = @(xdata, tdata) pd.ratefun(staticParams, alphafun(xdata{end}).*xdata{1});
end

for j = 1:length(btd)
    [sim(j).rate, sim(j).turnEti, sim(j).runEti] = simulateTurnTrain (btd(j), xnames, ratefun, pd.rateTx , nlarvae, 1);
    for k = 1:(length(xnames) - 1)
        ind = btd(j).findField(xnames{k});
        xv = btd(j).glt(ind).derivationMethod(pd.rateTx,  btd(j).glt(ind).xData,  btd(j).glt(ind).yData);
        sim(j).turnX(:,k) = interp1(pd.rateTx, xv, sim(j).turnEti);
        sim(j).runX(:,k) = interp1(pd.rateTx, xv, sim(j).runEti);
    end
    if (strcmpi(timeField, 'eti'))
        sim(j).turnT = sim(j).turnEti;
        sim(j).runT = sim(j).runEti;
    else
        ind = btd(j).findField([gq_switch '_' timeField]);
        tf =  btd(j).glt(ind).derivationMethod(pd.rateTx,  btd(j).glt(ind).xData,  btd(j).glt(ind).yData);
        sim(j).turnT = interp1(pd.rateTx, tf, sim(j).turnEti);
        sim(j).runT = interp1(pd.rateTx, tf, sim(j).runEti);
    end
    sim(j).turnExpnum = j * ones(size(sim(j).turnEti));
    sim(j).runExpnum = j * ones(size(sim(j).runEti));
end
pd.sim = sim;

pd.turnT = vertcat(sim.turnT);
pd.turnEti = vertcat(sim.turnEti);
pd.turnX = vertcat(sim.turnX);
pd.turnExpnum = vertcat(sim.turnExpnum);

pd.runT = vertcat(sim.runT);
pd.runEti = vertcat(sim.runEti);
pd.runX = vertcat(sim.runX);
pd.runExpnum = vertcat(sim.runExpnum);

existsAndDefault('deltaT', median(diff(opstruct.taxis)));
pd.deltaT = deltaT;

if (switchstim > 0) 
    period = opstruct.stim(switchstim).period;
else
    period = 0;
end
if (strcmpi(timeField, 'ton') || strcmpi(timeField, 'toff'))
   pd.tx = 0:pd.deltaT:period;   
   pd.pad = true;
   if (~exist('trange', 'var') || isempty('trange'))
       minTime = min(run(1).eti(run(1).ton >= 0 & run(1).toff >= 0));
       trange = minTime:period:max(run(1).eti + pd.deltaT); trange = trange([1 end]);
   end
    trange(1) = max(trange(1),min(pd.runEti));
    trange(2) = min(trange(2),max(pd.runEti));
else
    existsAndDefault('trange', [min(opstruct.taxis) max(opstruct.taxis)]);
    trange(1) = max(trange(1),min(pd.runEti));
    trange(2) = min(trange(2),max(pd.runEti));
    pd.tx = (trange(1):pd.deltaT:trange(end))';
    pd.pad = false;
end
existsAndDefault ('exprange', [min(pd.runExpnum) max(pd.runExpnum)]);
pd.exprange = exprange;

pd.timeRange = trange;
pd.etx = linspace(min(trange), max(trange), 100);

ti = pd.turnEti >= min(trange) & pd.turnEti < max(trange) & pd.turnExpnum >= min(exprange) & pd.turnExpnum <= max(exprange);
pd.turnT = pd.turnT(ti);
pd.turnEti = pd.turnEti(ti);
pd.turnX = pd.turnX(ti,:);
pd.turnExpnum = pd.turnExpnum(ti);

ri = pd.runEti >= min(trange) & pd.runEti < max(trange) & pd.runExpnum >= min(exprange) & pd.runExpnum <= max(exprange);
pd.runT = pd.runT(ri);
pd.runEti = pd.runEti(ri);
pd.runX = pd.runX(ri,:);
pd.runExpnum = pd.runExpnum(ri);

pd.params_0(end) = log(length(pd.turnT)/(length(pd.runT)*pd.deltaT));

pd.Q_alpha = 0.2*eye(ndim);
pd.Q_s = pd.Q_alpha;

pd.tparams_0 = [0 1000];
pd.alpha_0 = ones(length(pd.tx), ndim);
pd.sigma_0 = 10*ones(length(pd.tx), ndim);
pd.v0 = 1;
pd.w_0 = pd.Q_alpha*1000*eye(ndim);
pd.maxreps = 3;

pd.separateExperiments = false;
pd.period = period;
pd.tshift = 0; %tshift is incorrect for now -- fix later
pd.temporalratemod = [];
problemDescription = pd;

end

function hess = twobytwohess (x, xdata)
hess = zeros(size(xdata,1), 2, 2);
hess(:,1,1) = 2*x(1);
hess(:,1,2) = x(2);
hess(:,2,1) = x(2);
hess(:,2,2) = x(3);
end

