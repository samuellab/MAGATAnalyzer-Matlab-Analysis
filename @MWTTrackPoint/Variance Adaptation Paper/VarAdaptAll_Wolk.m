
function btdstruct = VarAdaptAll_Wolk(btdstruct, var, opstruct, Q_alpha, tau_alpha, deltaT, ppf, bayes, simtype, rescale_type, timeType, varargin)
%function btdstruct = VarAdaptAll(btdstruct, var, opstruct, Q_alpha, tau_alpha, deltaT, ppf, bayes, simtype, rescale_type, timeType)

% opstruct -- contains info about the experiment 

% var -- either [], or btdstruct.var
% Q_alpha -- variance of the prior for Point Process Filter
% tau_alpha -- correlation time for the Bayesian estimator (variance of the Gaussian prior = 2*deltaT/tau_alpha)
% deltaT -- time-bin width for PPF and Bayesian estimates of gain-rescaling factors
% ppf -- 0 or 1, runs PPF estimate of Alpha(t) or not
% bayes -- 0 or 1, runs Bayesian estimate of Alpha(t) or not
% simtype -- 'LowHigh' simulates the same number of experiments as are in btdstruct, with a rate function switching from low to high instantly
% simtype -- 'All' does the same, but uses the continuously changing rate function as calculated by Bayesian or PPF estimators
% simtype -- [] doesn't simulate anything
% rescale_type -- 'input' or 'output' ('output' might not work very well)
% timeType -- 'eti' calculates estimates vs. eti and then converts to cycle time
% timeType -- 'ton' calculates directly vs cycle time
% to include TTA bootstrap estimates add 'bootstrap', 1

bootstrap = false;
varargin = assignApplicable(varargin);

logdat = 1; % determines whether logP(data|model) or P(data|model) is used in Bayesian estimator -- shouldn't make a difference
nsim = 1; % how many times the data is simulated ( n=1 means the same number of experiments are simulated as were performed, with #larvae=mean(#larvae) )



if (isempty(var))
    % use prepVarianceSwitchingAnalysis if kernels need recalculating (not light & odor stimuli, e.g 20C08 or other neuron)
    btdstruct = BehaviorTriggeredData.prepVarianceSwitchingAnalysis_Wolk(btdstruct, opstruct);
    var = btdstruct.var;
    if(isfield(btdstruct, 'var_uv'))
        var_uv = btdstruct.var_uv;
    end
end

if(isfield(btdstruct, 'var_uv'))
    var_uv = btdstruct.var_uv;
end

pdegree = opstruct.pdegree;
if(size(var, 2)>=2)
    nstim = 2;
else
    nstim = 1;
end

if(~isfield(btdstruct, 'varops'))
    fn = {'stim', 'numLxBins', 'adaptationTime', 'Q_alpha', 'numTimeBins', 'timeRange', 'kernelTime', 'kernelDt'};
    for j = 1:length(fn)
        btdstruct.varops.(fn{j}) = opstruct.(fn{j});
    end
end

whichStim = find([btdstruct.varops.stim.switch]);

% TTA Kernels
if(~isfield(btdstruct, 'Kernels'))
    btdstruct = measureKernelsVsVar (btdstruct, var, bootstrap);
end

% Head Sweeps
% if(~isfield(btdstruct, 'HSavgs') && strcmpi(btdstruct.varops.stim(whichStim).ramptype, 'square'))
%     existsAndDefault('HSdt',1.5);
%     btdstruct =  measureHSstats (btdstruct, var, HSdt);
% end

% Rates vs. Stim.
if(~isfield(btdstruct, 'Rates'))
    btdstruct = measureRatesVsVar (btdstruct, var);
end
if(isfield(btdstruct, 'var_uv') && ~isfield(btdstruct, 'RatesU'))
    [~, btdstruct.RatesU] = measureRatesVsVar (btdstruct, var_uv);
end

% Rate vs Time
if(~isfield(btdstruct, 'RateVTime'))
    rateDt = 1;
    [~, btdstruct.RateVTime] = measureRateVtime(var(1), deltaT, rateDt);
end


% PPF & Bayes Estimates of Gain-Rescaling

if(bayes)
    Ddt_alpha = deltaT./tau_alpha;
    if(strcmpi(timeType, 'eti'))
        if(isfinite(btdstruct.varops.timeRange(2)))
            tx_bayes = btdstruct.varops.timeRange(1):deltaT:btdstruct.varops.timeRange(2);
        else
            tx_bayes = btdstruct.varops.timeRange(1):deltaT:max(btdstruct.var(1).fullensemble.eti);
        end
    else
        tx_bayes = 0:deltaT:btdstruct.var.period;
    end
end

if(ppf && ~isfield(btdstruct, 'alphavseti'))
    disp('starting PPF estimate of Data (Odor, Light)')
    btdstruct = runPPF(btdstruct, [], var, opstruct, Q_alpha, deltaT, timeType);
end

if( bayes && ~isfield(btdstruct, 'var_uv') && ~isfield(btdstruct, 'BayesAlphavsEti_OL'))
    tic
    if(nstim==2 && pdegree==2)
        disp('starting Bayesian estimate of Data (Odor, Light) - added rates') % r(xo,xl) = a*u^2 + b*u + c , with u=cos(theta)*xo+sin(theta)*xl
        btdstruct.BayesAlphavsEti_OL.add = BayesianModelEstimate([], btdstruct, var, opstruct, timeType, deltaT, Ddt_alpha, tx_bayes, logdat, [], rescale_type, 'add');
        disp('starting Bayesian estimate of Data (Odor, Light) - multiplied rates')% r(xo,xl) = r(xo)*r(xl), where r(x)=ax^2+b*x+c
        btdstruct.BayesAlphavsEti_OL.multiply = BayesianModelEstimate([], btdstruct, var, opstruct, timeType, deltaT, Ddt_alpha, tx_bayes, logdat, [], rescale_type, 'mult');
    else
        disp('starting Bayesian estimate of Data (Odor, Light)')
        btdstruct.BayesAlphavsEti_OL = BayesianModelEstimate([], btdstruct, var, opstruct, timeType, deltaT, Ddt_alpha, tx_bayes, logdat, [], rescale_type, []);
    end
    toc
end

if(isfield(btdstruct, 'var_uv') && ppf && ~isfield(btdstruct, 'alphavseti_UV'))
    tic
    disp('starting PPF estimate of Rotated Data - UV')
    [~, btdstruct.alphavseti_UV] = runPPF(btdstruct, [], var_uv, opstruct, Q_alpha, deltaT, timeType);
    toc
end
if(isfield(btdstruct, 'var_uv') && ppf && ~isfield(btdstruct, 'alphavseti_U'))
    tic
    disp('starting PPF estimate of Rotated Data - U')
    [~, btdstruct.alphavseti_U] = runPPF(btdstruct, [], var_uv(1), opstruct, Q_alpha, deltaT, timeType);
    toc
end


if( isfield(btdstruct, 'var_uv') && ~isempty(whichStim) && bayes )
    
    if(~isfield(btdstruct, 'BayesAlphavsEti_UV'))
        disp('starting Bayesian estimate of Rotated Data - UV')
        tic
        btdstruct.BayesAlphavsEti_UV = BayesianModelEstimate([], btdstruct, var_uv, opstruct, timeType, deltaT, Ddt_alpha, tx_bayes, logdat, [], rescale_type, []);
        toc
    end
    if(~isfield(btdstruct, 'BayesAlphavsEti_U'))
        disp('starting Bayesian estimate of Rotated Data - U')
        tic
        btdstruct.BayesAlphavsEti_U = BayesianModelEstimate([], btdstruct, var_uv(1), opstruct, timeType, deltaT, Ddt_alpha(1), tx_bayes, logdat, [], rescale_type, []);
        toc
    end
%     if(~isfield(btdstruct, 'BayesAlphavsEti_V'))
%         disp('starting Bayesian estimate of Rotated Data - V')
%         tic
%         btdstruct.BayesAlphavsEti_V = BayesianModelEstimate([], btdstruct, var_uv(2), opstruct, timeType, deltaT, Ddt_alpha(2), tx_bayes, logdat, [], rescale_type);
%         toc
%     end
end

% Scaled Rates
if( any([opstruct.stim.switch]) )
    if(  isfield(btdstruct, 'var_uv') && (bayes || ppf) && ~isfield(btdstruct, 'ScaledRatesU') )
        [~, btdstruct.ScaledRatesU] = adaptedRatesVsVar (btdstruct, btdstruct.var_uv(1).period, rescale_type);
    elseif( isfield(btdstruct, 'var') && (bayes || ppf) && ~isfield(btdstruct, 'ScaledRates'))
        [~, btdstruct.ScaledRates] = adaptedRatesVsVar (btdstruct, btdstruct.var(1).period, rescale_type);
    end
end

% DKL
if(~isfield(btdstruct, 'DKL'))
    [~, btdstruct.DKL] = measureDKLVsVar(btdstruct, var);
end
if(isfield(btdstruct, 'var_uv') && ~isfield(btdstruct, 'DKL_UV'))
    [~, btdstruct.DKL_UV] = measureDKLVsVar(btdstruct, var_uv);
end

% calculate fit params for alpha(sigma)=1/(sigma^2+sigma_0^2);
if(~isempty(whichStim) && nstim==1)
%     if(strcmpi(btdstruct.varops.stim(whichStim).ramptype, 'triangle'))
        if(ppf && ~isfield(btdstruct, 'gainrescaling_PPF'))
            btdstruct.gainrescaling_PPF = fitVarianceRescaling(btdstruct, opstruct, 'PPF', opstruct.stim(whichStim).ramptype, whichStim);
        end
        if(bayes && ~isfield(btdstruct, 'gainrescaling_Bayes'))
            btdstruct.gainrescaling_Bayes = fitVarianceRescaling(btdstruct, opstruct, 'Bayes', opstruct.stim(whichStim).ramptype, whichStim);
        end
        
%     end
end


% Simulated Data & Estimates

if(~isempty(whichStim))
    if(strcmpi(simtype, 'LowHigh') || strcmpi(simtype, 'All'))
        
        if(~isfield(btdstruct, 'SimData'))
            tic
            btdstruct.SimData = SimulateData(btdstruct, var, opstruct, 'ton', nsim, simtype); %keep 'ton' to have both runTon and runT in the simulated data
            toc
        end
        
        if(isfield(btdstruct, [simtype 'Sim']))
            eval(['Sim = btdstruct.' simtype 'Sim;']);
            
            if( ppf && ~isfield(Sim, 'alphavseti'))
                disp('starting PPF estimate of simulated data')
                for i=1:length(btdstruct.SimData)
                    [~, Sim.alphavseti{i}] = runPPF(btdstruct, btdstruct.SimData{i}, var, opstruct, Q_alpha, deltaT, timeType);
                    alpha_ton(i,:) = Sim.alphavseti{i}.alpha_ton;
                    valpha_ton(i,:) = Sim.alphavseti{i}.valpha_ton;
                end
                Sim.tx_ton_ppf = Sim.alphavseti{1}.tx_ton;
                Sim.alpha_ton_ppf = mean(alpha_ton);
                Sim.valpha_ton_ppf = mean(valpha_ton)./nsim;
            end
            
            if( bayes && ~isfield(Sim, 'BayesAlphavsEti'))
                disp('starting Bayesian estimate of simulated data')
                rP = btdstruct.BayesAlphavsEti.RateParams;
                for i=1:length(btdstruct.SimData)
                    tic
                    Sim.BayesAlphavsEti{i} = BayesianModelEstimate(btdstruct.SimData{i}, btdstruct, var, opstruct, timeType, deltaT, Ddt_alpha, tx_bayes, logdat, [], rescale_type, []);
                    alpha_ton(i,:) = Sim.BayesAlphavsEti{i}.alpha_ton;
                    valpha_ton(i,:) = Sim.BayesAlphavsEti{i}.valpha_ton;
                    toc
                end
                Sim.tx_ton_bayes = Sim.BayesAlphavsEti{1}.tx_ton;
                Sim.alpha_ton_bayes = mean(alpha_ton)/mean(mean(alpha_ton));
                Sim.valpha_ton_bayes = mean(valpha_ton)./nsim;
            end
            eval(['btdstruct.' simtype 'Sim = Sim;']);
        else
            if(ppf || bayes)
                if( ppf )
                    disp('starting PPF estimate of simulated data')
                    for i=1:length(btdstruct.SimData)
                        [~, Sim.alphavseti{i}] = runPPF(btdstruct, btdstruct.SimData{i}, var, opstruct, Q_alpha, deltaT, timeType);
                    end

                end
                if(bayes)
                    disp('starting Bayesian estimate of simulated data')
                    for i=1:length(btdstruct.SimData)
                        tic
                        Sim.BayesAlphavsEti{i} = BayesianModelEstimate(btdstruct.SimData{i}, btdstruct, var, opstruct, timeType, deltaT, Ddt_alpha, tx_bayes, logdat, [], rescale_type, []);
                        toc
                    end

                end
                eval(['btdstruct.' simtype 'Sim = Sim;']);
            end
        end
        
    end
end

end




function [btdstruct, DKL] = measureDKLVsVar (btdstruct, var)
opstruct = btdstruct.varops;
switchStim = find([opstruct.stim.switch]);

if(switchStim>0)
    
    nstim = size(var, 2);

    for j = 1:nstim
        
        deltaT = median(diff(var(j).fullensemble.eti));
        adata = var(j).noturn.x_conv;
        tdata = var(j).turn.x_conv;
        tton = var(j).turn.ton;
        ttoff = var(j).turn.toff;
        teti = var(j).turn.eti;
        rton = var(j).noturn.ton;
        rtoff = var(j).noturn.toff;
        reti = var(j).noturn.eti;
        
        if(strcmpi(opstruct.stim(switchStim).ramptype, 'square'))
            
            tH = tton > opstruct.adaptationTime & tton < ttoff & teti >= min(opstruct.timeRange) & teti < max(opstruct.timeRange);
            tL = ttoff > opstruct.adaptationTime & ttoff < tton & teti >= min(opstruct.timeRange) & teti < max(opstruct.timeRange);
            tU = tton > 0 & tton < opstruct.adaptationTime & teti >= min(opstruct.timeRange) & teti < max(opstruct.timeRange);
            tD = ttoff > 0  & ttoff < opstruct.adaptationTime & teti >= min(opstruct.timeRange) & teti < max(opstruct.timeRange);
            rH = rton > opstruct.adaptationTime & rton < rtoff & reti >= min(opstruct.timeRange) & reti < max(opstruct.timeRange);
            rL = rtoff > opstruct.adaptationTime & rtoff < rton & reti >= min(opstruct.timeRange) & reti < max(opstruct.timeRange);
            rU = rton > 0 & rton < opstruct.adaptationTime & reti >= min(opstruct.timeRange) & reti < max(opstruct.timeRange);
            rD = rtoff > 0  & rtoff < opstruct.adaptationTime & reti >= min(opstruct.timeRange) & reti < max(opstruct.timeRange);
        else
            tH = tton > opstruct.adaptationTime & tton > (3*max(tton)/8) & tton < (5*max(tton)/8) & teti >= min(opstruct.timeRange) & teti < max(opstruct.timeRange);
            tL = ttoff > opstruct.adaptationTime & ttoff > (3*max(ttoff)/8) & ttoff < (5*max(ttoff)/8) & teti >= min(opstruct.timeRange) & teti < max(opstruct.timeRange);
            tU = tton > opstruct.adaptationTime & tton < ttoff & teti >= min(opstruct.timeRange) & teti < max(opstruct.timeRange);
            tD = ttoff > opstruct.adaptationTime & ttoff < tton & teti >= min(opstruct.timeRange) & teti < max(opstruct.timeRange);
            rH = rton > opstruct.adaptationTime & rton > (3*max(rton)/8) & rton < (5*max(rton)/8) & reti >= min(opstruct.timeRange) & reti < max(opstruct.timeRange);
            rL = rtoff > opstruct.adaptationTime & rtoff > (3*max(rtoff)/8) & rtoff < (5*max(rtoff)/8) & reti >= min(opstruct.timeRange) & reti < max(opstruct.timeRange);
            rU = rton > opstruct.adaptationTime & rton < rtoff & reti >= min(opstruct.timeRange) & reti < max(opstruct.timeRange);
            rD = rtoff > opstruct.adaptationTime & rtoff < rton & reti >= min(opstruct.timeRange) & reti < max(opstruct.timeRange);
        end
        
        saH = std(adata(rH)); saL = std(adata(rL)); saU = std(adata(rU)); saD = std(adata(rD));
        stH = std(tdata(tH)); stL = std(tdata(tL)); stU = std(tdata(tU)); stD = std(tdata(tD));
        maH = mean(adata(rH)); maL = mean(adata(rL)); maU = mean(adata(rU)); maD = mean(adata(rD));
        mtH = mean(tdata(tH)); mtL = mean(tdata(tL)); mtU = mean(tdata(tU)); mtD = mean(tdata(tD));
        
        dkl1D(j).High = 0.5*(stH^2/saH^2 + (maH-mtH)^2/saH^2 - 1 + log(saH^2/stH^2));
        dkl1D(j).Low = 0.5*(stL^2/saL^2 + (maL-mtL)^2/saL^2 - 1 + log(saL^2/stL^2));
        dkl1D(j).Up = 0.5*(stU^2/saU^2 + (maU-mtU)^2/saU^2 - 1 + log(saU^2/stU^2));
        dkl1D(j).Down = 0.5*(stD^2/saD^2 + (maD-mtD)^2/saD^2 - 1 + log(saD^2/stD^2));
        
        
    end
    DKL.dkl1D = dkl1D;
    
    if(nstim==2)
        for k=1:2
            rdH(k,:) = var(k).noturn.x_conv(rH); rdL(k,:) = var(k).noturn.x_conv(rL);
            rdU(k,:) = var(k).noturn.x_conv(rU); rdD(k,:) = var(k).noturn.x_conv(rD);
            tdH(k,:) = var(k).turn.x_conv(tH); tdL(k,:) = var(k).turn.x_conv(tL);
            tdU(k,:) = var(k).turn.x_conv(tU); tdD(k,:) = var(k).turn.x_conv(tD);
        end

        covA = cov(rdH');
        covT = cov(tdH');
        muA = mean(rdH');
        muT = mean(tdH');
        a = covA\covT;
        b = covA\(muA-muT)';
        dkl2D.High = .5 * ( trace(a) + (muA-muT)*b - 2 + log(det(covA)/det(covT)) );
        
        covA = cov(rdL');
        covT = cov(tdL');
        muA = mean(rdL');
        muT = mean(tdL');
        a = covA\covT;
        b = covA\(muA-muT)';
        dkl2D.Low = .5 * ( trace(a) + (muA-muT)*b - 2 + log(det(covA)/det(covT)) );
        
        covA = cov(rdU');
        covT = cov(tdU');
        muA = mean(rdU');
        muT = mean(tdU');
        a = covA\covT;
        b = covA\(muA-muT)';
        dkl2D.Up = .5 * ( trace(a) + (muA-muT)*b - 2 + log(det(covA)/det(covT)) );
        
        covA = cov(rdD');
        covT = cov(tdD');
        muA = mean(rdD');
        muT = mean(tdD');
        a = covA\covT;
        b = covA\(muA-muT)';
        dkl2D.Down = .5 * ( trace(a) + (muA-muT)*b - 2 + log(det(covA)/det(covT)) );
        
        DKL.dkl2D = dkl2D;
    end
    
else
    
    
    m4 = @(x) mean(x)^4 + 6*mean(x)^2*std(x)^2 + 3*std(x)^4;
    sev = @(x) sqrt( (m4(x) - (length(x)-3)*std(x)^4/(length(x)-1)) / length(x) );
    for j = 1:size(var, 2)
        
        data = var(j);
        
        
        
        adata = data.noturn.x_conv;
        tdata = data.turn.x_conv;
        sa = std(adata);
        va = std(adata)^2;
        st = std(tdata);
        vt = std(tdata)^2;
        sva = sqrt(2*va^2/(length(adata)-1));
        svt = sqrt(2*vt^2/(length(tdata)-1));
        
%         sva = sev(adata);
%         svt = sev(tdata);
        
        ma = mean(adata);
        mt = mean(tdata);
        sma = sqrt(va/length(adata));
        smt = sqrt(vt/length(tdata));
        
        DKL.mean(j) = 0.5*(vt/va + (ma-mt)^2/va - 1 + log(va/vt));
        DKL.var(j) = .5*sva^2 * ( (va-vt-(mt-ma)^2)/va^2 ) + .5*svt^2 * (1/va - 1/vt) + sma^2 * (ma-mt)/va + smt^2 * (mt-ma)/va; 
        
        
    end
end
end

function SimulatedData = SimulateData(btdstruct, var, opstruct, timeField, nreps, simtype)

if(isfield(btdstruct, 'BayesAlphavsEti'))
    alphastruct = btdstruct.BayesAlphavsEti;
elseif(isfield(btdstruct, 'BayesAlphavsEti_OL'))
    alphastruct = btdstruct.BayesAlphavsEti_OL;
else
    alphastruct = btdstruct.alphavseti;
end

alphastruct.staticParams = alphastruct.params;

if(isfield(btdstruct, 'SimData'))
    SimulatedData = btdstruct.SimData;
    return
end

disp(['simulating ' num2str(nreps) ' set of experiments']);

deltaT = median(diff(var.fullensemble.eti));
pd = BehaviorTriggeredData.createProblemStructForRateFunFitting (var, timeField, opstruct.pdegree, opstruct.timeRange, [], deltaT);
pd.maxreps = 2;


for j = 1:length(btdstruct.btd), numan(j) = btdstruct.btd(j).es.numAnimals; end
Nlarvae = round(mean(numan));
% Nlarvae = 40;

dt = mean(diff(alphastruct.tx_ton));
nt = round(10/dt);

alpha_ton = alphastruct.alpha_ton/mean(alphastruct.alpha_ton);
alow = mean(alpha_ton(1+nt:round(end/2)));
ahigh = mean(alpha_ton(round(end/2)+nt:end));

if(strcmpi(simtype, 'LowHigh'))
    
    for j=1:nreps
        SimulatedData{j} = BehaviorTriggeredData.simulateVarSwitchFromPrevFit(btdstruct, alphastruct, pd, alow, ahigh, Nlarvae);
    end
    
elseif(strcmpi(simtype, 'All'))
    
    
    pd.tx_eti = alphastruct.tx;
    pd.alpha_eti = alphastruct.alpha;
    
    for k=1:nreps
        SimulatedData{k} = BehaviorTriggeredData.simulateVarSwitchFromPrevFit(btdstruct, alphastruct, pd, alow, ahigh, Nlarvae);
    end
    
end

end

function gainrescaling = fitVarianceRescaling(btdstruct, opstruct, estimator, type, switchstim)

rescalingfun = @(x,xdata) 1./sqrt(xdata.^2 + x^2) / mean(1./sqrt(xdata.^2 + x^2));
% rescalingfun = @(x,xdata) x(1)./sqrt(xdata.^2 + x(2)^2);

if(strcmpi(estimator, 'PPF'))
    if(length(btdstruct.varops.stim)>2)
        tx_ton = btdstruct.alphavseti_U.tx_ton;
        alpha_s = btdstruct.alphavseti_U.alpha_s_ton;
        valpha_s = squeeze(btdstruct.alphavseti_U.valpha_s_ton)';
    else
        tx_ton = btdstruct.alphavseti.tx_ton;
        alpha_s = btdstruct.alphavseti.alpha_s_ton(switchstim, :);
        valpha_s = squeeze(btdstruct.alphavseti.valpha_s_ton(switchstim, switchstim, :))';
    end
else
    if(length(btdstruct.varops.stim)>2)
        tx_ton = btdstruct.BayesAlphavsEti_U.tx_ton;
        alpha_s = btdstruct.BayesAlphavsEti_U.alpha_ton;
        valpha_s = btdstruct.BayesAlphavsEti_U.valpha_ton;
    else
        tx_ton = btdstruct.BayesAlphavsEti_OL.tx_ton;
        alpha_s = btdstruct.BayesAlphavsEti_OL.alpha_ton(switchstim, :);
        valpha_s = btdstruct.BayesAlphavsEti_OL.valpha_ton(switchstim, :);
    end
end

if(isfield(btdstruct, 'var_uv'))
    var = btdstruct.var_uv(1);
elseif(isfield(btdstruct, 'var'))
    var = btdstruct.var(1);
else
    disp('no var field to calculate signal variance');
    return
end

trange = opstruct.timeRange;
adaptTime = opstruct.adaptationTime;
teti = var.fullensemble.eti;
ton = var.fullensemble.ton;
xo = var.fullensemble.x_conv;

valid = teti>trange(1) & teti<trange(2);
adapted = (ton>adaptTime & ton<round(max(ton)/2)) | ton>(round(max(ton)/2)+adaptTime);

if(strcmpi(type, 'triangle'))
    [~,~,~,son_all] = meanyvsx(ton(valid), xo(valid), binEdgesFromCenters(tx_ton));
    % son_all = son_all';
    params.up = lsqcurvefit(rescalingfun, [1], son_all(tx_ton<(max(tx_ton)/2)), alpha_s(tx_ton<(max(tx_ton)/2)));
    params.down = lsqcurvefit(rescalingfun, [1], son_all(tx_ton>(max(tx_ton)/2)), alpha_s(tx_ton>(max(tx_ton)/2)));
    params.all = lsqcurvefit(rescalingfun, [1], son_all, alpha_s);
elseif(strcmpi(type, 'square'))
    
    high = tx_ton>=adaptTime & tx_ton<round(max(tx_ton)/2);
    low = tx_ton>(round(max(tx_ton)/2)+adaptTime);
    all = (tx_ton>adaptTime & tx_ton<round(max(tx_ton)/2)) | tx_ton>(round(max(tx_ton)/2)+adaptTime);
    
%     [~,~,~,son_all] = meanyvsx(ton(valid & adapted), xo(valid & adapted), binEdgesFromCenters(tx_ton));
%     params.high = lsqcurvefit(rescalingfun, [1 1], son_all(high), alpha_s(high));
%     params.low = lsqcurvefit(rescalingfun, [1 1], son_all(low), alpha_s(low));
    
    [~,~,~,son_all] = meanyvsx(ton(valid & adapted), xo(valid & adapted), binEdgesFromCenters(tx_ton));
    params.all = lsqcurvefit(rescalingfun, [1], son_all(all), alpha_s(all));
end
gainrescaling.estimator = estimator;
gainrescaling.fitfun = rescalingfun;
gainrescaling.son_all = son_all;
gainrescaling.params = params;
end


function [btdstruct, alphavseti, alphavston, scaledValues] = runPPF(btdstruct, pd, var, opstruct, Q_alpha, deltaT, timeType)


nstim = size(var, 2);

if(isempty(pd))
    pd = BehaviorTriggeredData.createProblemStructForRateFunFitting (var, 'eti', opstruct.pdegree, opstruct.timeRange, [], deltaT);
    pd.Q_alpha = Q_alpha * eye(nstim);
    pd.maxreps = 5;
end

pd.temporalratemod = [];


alphavseti = fitRateFunWithTemporalScalingND(pd);

%turn time estimate into time in cycle estimate

[alphavseti.alpha_ton, alphavseti.valpha_ton, alphavseti.tx_ton] = meanEstinCycle(alphavseti.alpha, alphavseti.valpha, alphavseti.tx, pd.tshift, pd.period, opstruct.timeRange, timeType);
[alphavseti.alpha_s_ton, alphavseti.valpha_s_ton, alphavseti.tx_ton] = meanEstinCycle(alphavseti.alpha_s, alphavseti.valpha_s, alphavseti.tx, pd.tshift, pd.period, opstruct.timeRange, timeType);

alphavseti.Q_alpha = Q_alpha;
btdstruct.alphavseti = alphavseti;


if(strcmpi(timeType, 'ton'))
    %directly calculate stretch factor vs cycle time
    pd = BehaviorTriggeredData.createProblemStructForRateFunFitting (var, 'ton', opstruct.pdegree, opstruct.timeRange, [], deltaT);
    
    if all(size(pd.alpha_0) == size(alphavseti.alpha_ton))
        pd.alpha_0 = alphavseti.alpha_ton;
    end
    pd.params_0 = alphavseti.staticParams;
    pd.tparams_0 = alphavseti.temporalParams;
    
    pd.Q_alpha = Q_alpha * eye(nstim) * pd.period / diff(pd.timeRange);
    [alphavston, scaledValues] = fitRateFunWithTemporalScalingND(pd);
    
    btdstruct.scaledValues = scaledValues;
    btdstruct.alphavston = alphavston;
end


end

function [btdstruct, ScaledRates] = adaptedRatesVsVar (btdstruct, period, rescale_type)
opstruct = btdstruct.varops;
switchStim = find([opstruct.stim.switch]);

if(isfield(btdstruct, 'BayesAlphavsEti_U'))
    scaledValues = btdstruct.BayesAlphavsEti_U.scaledValues;
elseif(isfield(btdstruct, 'BayesAlphavsEti_OL'))
    scaledValues = btdstruct.BayesAlphavsEti_OL.scaledValues;
elseif(isfield(btdstruct, 'BayesAlphavsEti'))
    scaledValues = btdstruct.BayesAlphavsEti.scaledValues;
else
    scaledValues = btdstruct.scaledValues;
end
tton = scaledValues.turnTon;
teti = scaledValues.turnEti;
rton = scaledValues.runTon;
reti = scaledValues.runEti;

tval = scaledValues.tval;
rval = scaledValues.rval;
deltaT = scaledValues.deltaT;

tvalid = teti >= min(opstruct.timeRange) & teti < max(opstruct.timeRange);
rvalid = reti >= min(opstruct.timeRange) & reti < max(opstruct.timeRange);

if(strcmpi(opstruct.stim(switchStim).ramptype, 'square'))
    thigh = tton > opstruct.adaptationTime & tton < period/2 & tvalid;
    tlow = tton > period/2 + opstruct.adaptationTime & tton < period & tvalid;
    tup = tton > 0 & tton < opstruct.adaptationTime & tvalid;
    tdown = tton > period/2  & tton < period/2 + opstruct.adaptationTime & tvalid;
    
    rhigh = rton > opstruct.adaptationTime & rton < period/2 & rvalid;
    rlow = rton > period/2 + opstruct.adaptationTime & rton < period & rvalid;
    rup = rton > 0 & rton < opstruct.adaptationTime & rvalid;
    rdown = rton > period/2  & rton < period/2 + opstruct.adaptationTime & rvalid;
else

    thigh = tton > opstruct.adaptationTime & tton < (period/2) & tvalid;
    tlow = tton > (period/2 + opstruct.adaptationTime) & tton < period & tvalid;
    tup = tton < opstruct.adaptationTime & tvalid;
    tdown = tton >(period/2) & tton<(period/2 + opstruct.adaptationTime) & tvalid;
    
    rhigh = rton > opstruct.adaptationTime & rton < (period/2) & rvalid;
    rlow = rton > (period/2 + opstruct.adaptationTime) & rton < period & rvalid;
    rup = rton < opstruct.adaptationTime & rvalid;
    rdown = rton >(period/2) & rton<(period/2 + opstruct.adaptationTime) & rvalid;
end

if(strcmpi(rescale_type, 'input'))
    for j = 1:size(rval, 2)
        lx = linspace(-3*std(rval(rvalid,j)), 3*std(rval(rvalid,j)), opstruct.numLxBins);
        ScaledRates.scaledRFHigh(j) = turnRates(tval(thigh, j), rval(rhigh, j), deltaT, [], opstruct);
        ScaledRates.scaledRFLow(j) = turnRates(tval(tlow, j), rval(rlow, j), deltaT, [], opstruct);
        ScaledRates.scaledRFUp(j) = turnRates(tval(tup, j), rval(rup, j), deltaT, [], opstruct);
        ScaledRates.scaledRFDown(j) = turnRates(tval(tdown, j), rval(rdown, j), deltaT, [], opstruct);
    end
    
elseif(strcmpi(rescale_type, 'output'))
    
    alphaT = scaledValues.alphaT;
    turnX = scaledValues.turnX;
    runX = scaledValues.runX;

    for j = 1:size(rval, 2)
        lx = linspace(-3*std(runX(rvalid,j)), 3*std(runX(rvalid,j)), opstruct.numLxBins);
        ScaledRates.scaledRFHigh(j) = scaledturnRates(turnX(thigh, j), alphaT(thigh, j), runX(rhigh, j), deltaT, [], opstruct);
        ScaledRates.scaledRFLow(j) = scaledturnRates(turnX(tlow, j), alphaT(tlow, j), runX(rlow, j), deltaT, [], opstruct);
        ScaledRates.scaledRFUp(j) = scaledturnRates(turnX(tup,j), alphaT(tup,j), runX(rup, j), deltaT, [], opstruct);
        ScaledRates.scaledRFDown(j) = scaledturnRates(turnX(tdown,j), alphaT(tdown,j), runX(rdown, j), deltaT, [], opstruct);
    end
end

btdstruct.scaledRFUp = ScaledRates.scaledRFUp;
btdstruct.scaledRFDown = ScaledRates.scaledRFDown;
btdstruct.scaledRFLow = ScaledRates.scaledRFLow;
btdstruct.scaledRFHigh = ScaledRates.scaledRFHigh;
end

function rf = turnRates2D(xt, xr, deltaT, tinds, rinds, nBins)

td(:,1) = xt(tinds, 1);
td(:,2) = xt(tinds, 2);
rd(:,1) = xr(rinds, 1);
rd(:,2) = xr(rinds, 2);
lx(:,1) = linspace(-3*std(rd(:,1)), 3*std(rd(:,1)), nBins);
lx(:,2) = linspace(-3*std(rd(:,2)), 3*std(rd(:,2)), nBins);
[lxx, lxy] = meshgrid(lx(:,1), lx(:,2));
lxdata = [lxx(:) lxy(:)];



h = makeIm(rd(:,1), rd(:,2), binEdgesFromCenters(lx(:,1)), binEdgesFromCenters(lx(:,2)));
ht = makeIm(td(:,1), td(:,2), binEdgesFromCenters(lx(:,1)), binEdgesFromCenters(lx(:,2)));
rhist = reshape(h, size(lxx));
thist = reshape(ht, size(lxx));

rf.lx = lx;
rf.rate(:,1) = 1/deltaT * sum(thist, 1)./sum(rhist, 1);
rf.rate_eb(:,1) = 1/deltaT * sqrt(sum(thist, 1))./sum(rhist, 1);
rf.rate(:,2) = 1/deltaT * sum(thist, 2)./sum(rhist, 2);
rf.rate_eb(:,2) = 1/deltaT * sqrt(sum(thist, 2))./sum(rhist, 2);


end
   
function rf = scaledturnRates (xt, alphaT, xr, deltaT, lx, opstruct)

op = optimoptions('fminunc');
op.Display = 'off';
op.Algorithm = 'quasi-newton';

if(isempty(lx))
    lx = linspace(-3*std(xr), 3*std(xr), opstruct.numLxBins);
end
existsAndDefault('lx', linspace(-3*std(xr), 3*std(xr), opstruct.numLxBins));

[nt, ~, tbins] = histcounts(xt, binEdgesFromCenters(lx));
nr = histcounts(xr, binEdgesFromCenters(lx));

for i=1:length(lx), at(i) = mean(alphaT(tbins==i)); sat(i) = std(alphaT(tbins==i)); end

rf.lx = lx;
rf.rate = 1/deltaT * at .* nt./nr;
rf.rate_eb = 1/deltaT * at .* sqrt(nt) ./ nr;

nlogp = @(p) -sum(polyval(p, xt)) + deltaT*sum(exp(polyval(p,xr)));
p0 = [mean(xt)/var(xr) log(length(xt)/(length(xr)*deltaT))];
rf.rateFitLin = fminunc(nlogp, p0, op);
p0 = [0 rf.rateFitLin];
rf.rateFitQuad = fminunc(nlogp, p0, op);

end 

function rf = turnRates (xt, xr, deltaT, lx, opstruct)

op = optimoptions('fminunc');
op.Display = 'off';
op.Algorithm = 'quasi-newton';

if(isempty(lx))
    lx = linspace(-3*std(xr), 3*std(xr), opstruct.numLxBins);
end
existsAndDefault('lx', linspace(-3*std(xr), 3*std(xr), opstruct.numLxBins));

nt = histcounts(xt, binEdgesFromCenters(lx));
nr = histcounts(xr, binEdgesFromCenters(lx));

rf.lx = lx;
rf.rate = 1/deltaT * nt ./ nr;
rf.rate_eb = 1/deltaT * sqrt(nt) ./ nr;

nlogp = @(p) -sum(polyval(p, xt)) + deltaT*sum(exp(polyval(p,xr)));
nlogp2 = @(p) -sum(polyval(p, xt)) -sum(log(deltaT*ones(1,length(xt)))) + deltaT*sum(exp(polyval(p,xr)));
p0 = [mean(xt)/var(xr) log(length(xt)/(length(xr)*deltaT))];
%p0 = [0.5 -1.6];
[rf.rateFitLin, nllh_lin] = fminunc(nlogp, p0, op);
p0 = [0 rf.rateFitLin];
[rf.rateFitQuad, nllh_quad] = fminunc(nlogp, p0, op);
rf.rateFitLin_llh = - nlogp2(rf.rateFitLin);
rf.rateFitQuad_llh = - nlogp2(rf.rateFitQuad);

end    

function [btdstruct, HS_stats] = measureHSstats (btdstruct, var, dt)

opstruct = btdstruct.varops;
switchStim = find([opstruct.stim.switch]);

taxis = -(-opstruct.kernelTime:0.5:opstruct.kernelTime);

if(switchStim>0)
    
    for j = 1:size(var, 2)
        maxT = 1.25;
        period = opstruct.stim(j).period;
        tBins = 0:dt:round(period/2-maxT);
        nt = dt/abs(median(diff(taxis)));
        
        firsths.name = 'hsnum';
        firsths.validop = @(x) x == 1;
        
        onehs.name = 'numhs';
        onehs.validop = @(x) x==1;
        
        tvalid.name = 'start_eti';
        tvalid.validop = @(x) x>= min(opstruct.timeRange) & x< max(opstruct.timeRange);
        
        lowvartime.name = [opstruct.stim(j).gqname '_var_toff'];
        
        
        acchs = [btdstruct.btd.acchs]';
        rejhs = [btdstruct.btd.rejhs]';
        acchseti = [acchs.start_eti]';
        rejhseti = [rejhs.start_eti]';
        
        accnhs = [acchs.hsnum]';
        rejnhs = [rejhs.hsnum]';
        firstacchs = accnhs == 1;
        firstrejhs = rejnhs == 1;
       

        delta_t_stim = median(diff(btdstruct.btd(1).all.eti));
        for i=1:(length(tBins)-1)
            
            lowvartime.validop = @(x) x > tBins(i) & x < tBins(i+1);

            newtaxis = 0:delta_t_stim:maxT;
            dm_turn = btdstruct.btd.behaviorTriggeredDataMatrix('turn', 'start', opstruct.stim(j).gqname, newtaxis, 'conditions', [onehs lowvartime tvalid]);
            dm_acchs = btdstruct.btd.behaviorTriggeredDataMatrix('acchs', 'start', opstruct.stim(j).gqname, newtaxis, 'conditions', [firsths lowvartime tvalid] );
            dm_rejhs = btdstruct.btd.behaviorTriggeredDataMatrix('rejhs', 'start', opstruct.stim(j).gqname, newtaxis, 'conditions', [firsths lowvartime tvalid]);
            
            meanDeltaAll = sum(dm_turn,2)*delta_t_stim/maxT;
            meanDeltaAll = meanDeltaAll(isfinite(meanDeltaAll)); 
            HS_stats.AllAvg(i) = mean(meanDeltaAll);
            HS_stats.AllAvg_err(i) = std(meanDeltaAll)/sqrt(length(meanDeltaAll));
            
            meanDeltaAcchs = sum(dm_acchs,2)*delta_t_stim/maxT;
            meanDeltaAcchs = meanDeltaAcchs(isfinite(meanDeltaAcchs)); 
            HS_stats.AccAvg(i) = mean(meanDeltaAcchs);
            HS_stats.AccAvg_err(i) = std(meanDeltaAcchs)/sqrt(length(meanDeltaAcchs)); %makes false assumption about independence of HS; ok if only 1 (or 2) first headsweeps per 1.5s of experiment
            
            meanDeltaRejhs = sum(dm_rejhs,2)*delta_t_stim/maxT;
            meanDeltaRejhs = meanDeltaRejhs(isfinite(meanDeltaRejhs)); 
            HS_stats.RejAvg(i) = mean(meanDeltaRejhs);
            HS_stats.RejAvg_err(i) = std(meanDeltaRejhs)/sqrt(length(meanDeltaRejhs));
            
            HS_stats.numAcc(i) = length(meanDeltaAcchs);
            HS_stats.numRej(i) = length(meanDeltaRejhs);
            
%             nhs = sum(isfinite(dm_acchs),1);
%             dm_acchs(~isfinite(dm_acchs)) = 0;
%             tta_acchs = sum(dm_acchs, 1)./nhs;
%             ttv_acchs =  sum( (dm_acchs-tta_acchs).^2, 1)./nhs;
%             
%             dm_rejhs(~isfinite(dm_rejhs)) = 0;
%             tta_rejhs = sum(dm_rejhs, 1)./sum(isfinite(dm_rejhs), 1);
%             ttv_rejhs = sum( (dm_rejhs-tta_rejhs).^2, 1 )./(sum(isfinite(dm_rejhs), 1)-1);
%             
%             HS_stats.tx = taxis;
%             HS_stats.AccTTA{i} = tta_acchs;
%             HS_stats.RejTTA{i} = tta_rejhs;
%             
%             HS_stats.SDall(i) = sqrt( sum(vt(taxis > 0 & taxis < 1.5))/(nt^2) );
%             
%             
%             HS_stats.AccAvg(i) = sum(tta_acchs(taxis > 0 & taxis < 1.5))/nt;
%             HS_stats.AccAvg_err(i) = sqrt( sum(ttv_acchs(taxis > 0 & taxis < 1.5))/(nt^2) );
%             HS_stats.RejAvg(i) = sum(tta_rejhs(taxis > 0 & taxis < 1.5));
%             HS_stats.RejAvg_err(i) = sqrt( sum(ttv_rejhs(taxis > 0 & taxis < 1.5))/(nt^2) );
%             
%             HS_stats.numAcc(i) = size(dm_acchs, 1);
%             HS_stats.numRej(i) = size(dm_rejhs, 1);
            
        end
        HS_stats.tt = tBins(1:end-1);
        btdstruct.HSavgs = HS_stats;
        
        
    end
end
end

function [btdstruct, Kernels] = measureKernelsVsVar (btdstruct, var, bootstrap)



opstruct = btdstruct.varops;
switchStim = find([opstruct.stim.switch]);

firsths.name = 'hsnum';
firsths.validop = @(x) x == 1;

if(switchStim>0)
    
    for j = 1:size(var, 2)
        
        period = var(j).period;
        
        tton = var(j).turn.ton;
        ttoff = var(j).turn.toff;
        teti = var(j).turn.eti;
        
        accton = var(j).acchs.ton;
        acctoff = var(j).acchs.toff;
        rejton = var(j).rejhs.ton;
        rejtoff = var(j).rejhs.toff;
        acceti = var(j).acchs.eti;
        rejeti = var(j).rejhs.eti;
        
        tt = [btdstruct.btd.turn]';
        turneti = [tt.start_eti]';
        turntnum = [tt.tnum]';
        acchs = [btdstruct.btd.acchs]';
        rejhs = [btdstruct.btd.rejhs]';
        acchseti = [acchs.start_eti]';
        rejhseti = [rejhs.start_eti]';
        
        accnhs = [acchs.hsnum]';
        rejnhs = [rejhs.hsnum]';
        firstacchs = accnhs == 1;
        firstrejhs = rejnhs == 1;

        turnvalid = turneti >= min(btdstruct.varops.timeRange) & turneti < max(btdstruct.varops.timeRange);
        acchsvalid = acchseti >= min(btdstruct.varops.timeRange) & acchseti < max(btdstruct.varops.timeRange);
        rejhsvalid = rejhseti >= min(btdstruct.varops.timeRange) & rejhseti < max(btdstruct.varops.timeRange);

        
        
        taxis = -(-opstruct.kernelTime:opstruct.kernelDt:opstruct.kernelTime);
        
        if strcmp(opstruct.stim(j).gqname, 'led1ValDiff')
            gqname_levels = ['led1Val'];
        elseif strcmp(opstruct.stim(j).gqname, 'led2ValDiff')
            gqname_levels = ['led2Val'];
        end
        
        if strcmp(opstruct.stim(j).gqname, 'led1ValDiff') || strcmp(opstruct.stim(j).gqname, 'led2ValDiff')
            dm_val_lvl = btdstruct.btd.behaviorTriggeredDataMatrix('turn', 'start', gqname_levels, taxis);
            dm_val_lvl = dm_val_lvl(turnvalid, :);
        end
        
        [dm_val,~,dm_expnum] = btdstruct.btd.behaviorTriggeredDataMatrix('turn', 'start', opstruct.stim(j).gqname, taxis);
        dm_val = dm_val(turnvalid, :);
        dm_expnum = dm_expnum(turnvalid, :);
        turntnum = turntnum(turnvalid, :);
        
        dm_acchs = btdstruct.btd.behaviorTriggeredDataMatrix('acchs', 'end', opstruct.stim(j).gqname, taxis);
        dm_acchs = dm_acchs(acchsvalid & firstacchs, :);
        
        dm_rejhs = btdstruct.btd.behaviorTriggeredDataMatrix('rejhs', 'end', opstruct.stim(j).gqname, taxis);
        dm_rejhs = dm_rejhs(rejhsvalid & firstrejhs, :);

                
        if(strcmpi(opstruct.stim(switchStim).ramptype, 'square'))
           
            tH = tton > opstruct.adaptationTime & tton < ttoff & teti >= min(opstruct.timeRange) & teti < max(opstruct.timeRange);
            tL = ttoff > opstruct.adaptationTime & ttoff < tton & teti >= min(opstruct.timeRange) & teti < max(opstruct.timeRange);
            tU = tton > 0 & tton < (period/2-opstruct.adaptationTime) & teti >= min(opstruct.timeRange) & teti < max(opstruct.timeRange);
            tD = ttoff > 0  & ttoff < (period/2-opstruct.adaptationTime) & teti >= min(opstruct.timeRange) & teti < max(opstruct.timeRange);
            
            accH = accton > opstruct.adaptationTime & accton < acctoff & acceti >= min(opstruct.timeRange) & acceti < max(opstruct.timeRange);
            accL = acctoff > opstruct.adaptationTime & acctoff < accton & acceti >= min(opstruct.timeRange) & acceti < max(opstruct.timeRange);
            accU = accton > 0 & accton < (period/2-opstruct.adaptationTime) & acceti >= min(opstruct.timeRange) & acceti < max(opstruct.timeRange);
            accD = acctoff > 0  & acctoff < (period/2-opstruct.adaptationTime) & acceti >= min(opstruct.timeRange) & acceti < max(opstruct.timeRange);
            rejH = rejton > opstruct.adaptationTime & rejton < rejtoff & rejeti >= min(opstruct.timeRange) & rejeti < max(opstruct.timeRange);
            rejL = rejtoff > opstruct.adaptationTime & rejtoff < rejton & rejeti >= min(opstruct.timeRange) & rejeti < max(opstruct.timeRange);
            rejU = rejton > 0 & rejton < (period/2-opstruct.adaptationTime) & rejeti >= min(opstruct.timeRange) & rejeti < max(opstruct.timeRange);
            rejD = rejtoff > 0  & rejtoff < (period/2-opstruct.adaptationTime) & rejeti >= min(opstruct.timeRange) & rejeti < max(opstruct.timeRange);
        else
            tH = tton > opstruct.adaptationTime & tton > (3*max(tton)/8) & tton < (5*max(tton)/8) & teti >= min(opstruct.timeRange) & teti < max(opstruct.timeRange);
            tL = ttoff > opstruct.adaptationTime & ttoff > (3*max(ttoff)/8) & ttoff < (5*max(ttoff)/8) & teti >= min(opstruct.timeRange) & teti < max(opstruct.timeRange);
            tU = tton > opstruct.adaptationTime & tton < ttoff & teti >= min(opstruct.timeRange) & teti < max(opstruct.timeRange);
            tD = ttoff > opstruct.adaptationTime & ttoff < tton & teti >= min(opstruct.timeRange) & teti < max(opstruct.timeRange);
        end
        
        if strcmp(opstruct.stim(j).gqname, 'led1ValDiff') || strcmp(opstruct.stim(j).gqname, 'led2ValDiff')
            dm = dm_val_lvl(tH, :);
            dm(~isfinite(dm)) = 0;
            Kernels.High_lvl(j).taxis = taxis;
            Kernels.High_lvl(j).tta = sum(dm)./sum(isfinite(dm_val(tH,:)));
            Kernels.High_lvl(j).tta = Kernels.High_lvl(j).tta - mean(btdstruct.btd.behaviorTriggeredAverage('turn','start', gqname_levels , 10:0.1:15));
            Kernels.High_lvl(j).nturns = nnz(tH);
        
            dm = dm_val_lvl(tL, :);
            dm(~isfinite(dm)) = 0;
            Kernels.Low_lvl(j).taxis = taxis;
            Kernels.Low_lvl(j).tta = sum(dm)./sum(isfinite(dm_val(tL,:)));
            Kernels.Low_lvl(j).tta = Kernels.Low_lvl(j).tta - mean(btdstruct.btd.behaviorTriggeredAverage('turn','start', gqname_levels , 10:0.1:15));
            Kernels.Low_lvl(j).nturns = nnz(tL);
        end
        
        dm = dm_val(tH, :);
        dm_e = dm_expnum(tH, :);
        dm_tnum = turntnum(tH, :);
        dm(~isfinite(dm)) = 0;
        Kernels.High(j).taxis = taxis;
        Kernels.High(j).tta = sum(dm)./sum(isfinite(dm_val(tH,:))); 
        Kernels.High(j).nturns = nnz(tH);
        Kernels.High(j).dm = dm;
        Kernels.High(j).dm_e = dm_e;
        Kernels.High(j).dm_tnum = dm_tnum;
% code for bootstrapping TTA randomly sampling from experiments and tracks:
        if (bootstrap == 1)
            [meanTTA, meanBS, stdBS] = bootstrapTTAfromBTDM(dm, dm_e, dm_tnum, taxis, 100);
            Kernels.High(j).bootstrapKer = meanTTA;   
            Kernels.High(j).bootstrapMean = meanBS;
            Kernels.High(j).bootstrapStd = stdBS;    
        end
        
        dm = dm_val(tL, :);
        dm_e = dm_expnum(tL, :);
        dm_tnum = turntnum(tL, :);
        dm(~isfinite(dm)) = 0;
        Kernels.Low(j).taxis = taxis;
        Kernels.Low(j).tta = sum(dm)./sum(isfinite(dm_val(tL,:)));
        Kernels.Low(j).nturns = nnz(tL);
        Kernels.Low(j).dm = dm;
        Kernels.Low(j).dm_e = dm_e;
        Kernels.Low(j).dm_tnum = dm_tnum;
% code for bootstrapping TTA randomly sampling from experiments and tracks:
        if (bootstrap == 1)
            [meanTTA, meanBS, stdBS] = bootstrapTTAfromBTDM(dm, dm_e, dm_tnum, taxis, 100);
            Kernels.Low(j).bootstrapKer = meanTTA;
            Kernels.Low(j).bootstrapMean = meanBS;
            Kernels.Low(j).bootstrapStd = stdBS;
        end
        
        dm = dm_val(tU, :);
        dm(~isfinite(dm)) = 0;
        Kernels.Up(j).taxis = taxis;
        Kernels.Up(j).tta = sum(dm)./sum(isfinite(dm_val(tU,:)));
        Kernels.Up(j).nturns = nnz(tU);
        
        dm = dm_val(tD, :);
        dm(~isfinite(dm)) = 0;
        Kernels.Down(j).taxis = taxis;
        Kernels.Down(j).tta = sum(dm)./sum(isfinite(dm_val(tD,:)));
        Kernels.Down(j).nturns = nnz(tD);
        
        if(strcmpi(opstruct.stim(switchStim).ramptype, 'square'))
            dm_acc = dm_acchs(accH, :);
            dm_acc(~isfinite(dm_acc)) = 0;
            dm_rej = dm_rejhs(rejH, :);
            dm_rej(~isfinite(dm_rej)) = 0;
            Kernels.High(j).tta_acchs = sum(dm_acc)./sum(isfinite(dm_acchs(accH, :)));
            Kernels.High(j).tta_rejhs = sum(dm_rej)./sum(isfinite(dm_rejhs(rejH, :)));
            dm_acc = dm_acchs(accL, :);
            dm_acc(~isfinite(dm_acc)) = 0;
            dm_rej = dm_rejhs(rejL, :);
            dm_rej(~isfinite(dm_rej)) = 0;
            Kernels.Low(j).tta_acchs = sum(dm_acc)./sum(isfinite(dm_acchs(accL, :)));
            Kernels.Low(j).tta_rejhs = sum(dm_rej)./sum(isfinite(dm_rejhs(rejL, :)));
            dm_acc = dm_acchs(accU, :);
            dm_acc(~isfinite(dm_acc)) = 0;
            dm_rej = dm_rejhs(rejU, :);
            dm_rej(~isfinite(dm_rej)) = 0;
            Kernels.Up(j).tta_acchs = sum(dm_acc)./sum(isfinite(dm_acchs(accU, :)));
            Kernels.Up(j).tta_rejhs = sum(dm_rej)./sum(isfinite(dm_rejhs(rejU, :)));
            dm_acc = dm_acchs(accD, :);
            dm_acc(~isfinite(dm_acc)) = 0;
            dm_rej = dm_rejhs(rejD, :);
            dm_rej(~isfinite(dm_rej)) = 0;
            Kernels.Down(j).tta_acchs = sum(dm_acc)./sum(isfinite(dm_acchs(accD, :)));
            Kernels.Down(j).tta_rejhs = sum(dm_rej)./sum(isfinite(dm_rejhs(rejD, :)));
        end
        
    end

else
    
    for j = 1:size(var, 2)
        teti = var(j).turn.eti;
        turnvalid = teti >= min(opstruct.timeRange) & teti < max(opstruct.timeRange);
        
        taxis = -(-2:opstruct.kernelDt:opstruct.kernelTime);
        dm_val = btdstruct.btd.behaviorTriggeredDataMatrix('turn', 'start', opstruct.stim(j).gqname, taxis);
        dm_conv = btdstruct.btd.behaviorTriggeredDataMatrix('turn', 'start', [opstruct.stim(j).gqname 'TurnLin'], 0);
        dm_val = dm_val(turnvalid, :);
        dm_val(~isfinite(dm_val)) = 0;
        dm_conv = dm_conv(turnvalid);
        Kernels(j).taxis = taxis;
        Kernels(j).tta = sum(dm_val)./sum(isfinite(dm_conv));
    end
    
end


btdstruct.Kernels = Kernels;


end

function [btdstruct, Rates] = measureRatesVsVar (btdstruct, var)
    

opstruct = btdstruct.varops;
switchStim = find([opstruct.stim.switch]);

op = optimoptions('fminunc');
op.Display = 'off';
op.Algorithm = 'quasi-newton';


if(switchStim>0)
     
    for j = 1:size(var, 2)

        period = var(j).period;
        
        deltaT = median(diff(var(j).fullensemble.eti));
        tton = var(j).turn.ton;
        ttoff = var(j).turn.toff;
        teti = var(j).turn.eti;
        rton = var(j).noturn.ton;
        rtoff = var(j).noturn.toff;
        reti = var(j).noturn.eti;
        
        if(strcmpi(opstruct.stim(switchStim).ramptype, 'square'))
           
            tH = tton > opstruct.adaptationTime & tton < ttoff & teti >= min(opstruct.timeRange) & teti < max(opstruct.timeRange);
            tL = ttoff > opstruct.adaptationTime & ttoff < tton & teti >= min(opstruct.timeRange) & teti < max(opstruct.timeRange);
            tU = tton > 0 & tton < opstruct.adaptationTime & teti >= min(opstruct.timeRange) & teti < max(opstruct.timeRange);
            tD = ttoff > 0  & ttoff < opstruct.adaptationTime & teti >= min(opstruct.timeRange) & teti < max(opstruct.timeRange); 
            rH = rton > opstruct.adaptationTime & rton < rtoff & reti >= min(opstruct.timeRange) & reti < max(opstruct.timeRange);
            rL = rtoff > opstruct.adaptationTime & rtoff < rton & reti >= min(opstruct.timeRange) & reti < max(opstruct.timeRange);
            rU = rton > 0 & rton < opstruct.adaptationTime & reti >= min(opstruct.timeRange) & reti < max(opstruct.timeRange);
            rD = rtoff > 0  & rtoff < opstruct.adaptationTime & reti >= min(opstruct.timeRange) & reti < max(opstruct.timeRange);
        else
            
            tH = tton>(period/2-opstruct.adaptationTime) & tton<(period/2+opstruct.adaptationTime) & teti >= min(opstruct.timeRange) & teti < max(opstruct.timeRange);
            tL = ttoff>(period/2-opstruct.adaptationTime) & ttoff<(period/2+opstruct.adaptationTime) & teti >= min(opstruct.timeRange) & teti < max(opstruct.timeRange);
            tU = tton > opstruct.adaptationTime & tton < ttoff & teti >= min(opstruct.timeRange) & teti < max(opstruct.timeRange);
            tD = ttoff > opstruct.adaptationTime & ttoff < tton & teti >= min(opstruct.timeRange) & teti < max(opstruct.timeRange);
            rH = rton>(period/2-opstruct.adaptationTime) & rton<(period/2+opstruct.adaptationTime) & reti >= min(opstruct.timeRange) & reti < max(opstruct.timeRange);
            rL = rtoff>(period/2-opstruct.adaptationTime) & rtoff<(period/2+opstruct.adaptationTime) & reti >= min(opstruct.timeRange) & reti < max(opstruct.timeRange);
            rU = rton > opstruct.adaptationTime & rton < rtoff & reti >= min(opstruct.timeRange) & reti < max(opstruct.timeRange);
            rD = rtoff > opstruct.adaptationTime & rtoff < rton & reti >= min(opstruct.timeRange) & reti < max(opstruct.timeRange);
            %rH = rton > opstruct.adaptationTime & rton > (7*max(rton)/16) & rton < (9*max(rton)/16) & reti >= min(opstruct.timeRange) & reti < max(opstruct.timeRange);
            %rL = rtoff > opstruct.adaptationTime & rtoff > (7*max(rtoff)/16) & rtoff < (9*max(rtoff)/16) & reti >= min(opstruct.timeRange) & reti < max(opstruct.timeRange);
        end
        
        xt = var(j).turn.x_conv(tH);
        xr = var(j).noturn.x_conv(rH);
        lx = linspace(-3*std(xr), 3*std(xr), opstruct.numLxBins);   
        Rates.rfHighVar(j) = turnRates (xt, xr, deltaT, lx, opstruct); 
             
        xt = var(j).turn.x_conv_HK(tH);
        xr = var(j).noturn.x_conv_HK(rH);
        lx = linspace(-3*std(xr), 3*std(xr), opstruct.numLxBins);   
        Rates.rfHighVar_HK(j) = turnRates (xt, xr, deltaT, lx, opstruct);
        
        if strcmp(opstruct.stim(j).gqname, 'led1ValDiff') || strcmp(opstruct.stim(j).gqname, 'led2ValDiff')
            xt = var(j).turn.x_conv_lvl(tH);
            xr = var(j).noturn.x_conv_lvl(rH);
            lx = linspace(mean(xr)-3*std(xr), mean(xr)+3*std(xr), opstruct.numLxBins);   
            Rates.rfHighVar_lvl(j) = turnRates (xt, xr, deltaT, lx, opstruct);
        
            xt = var(j).turn.x_conv_lvl(tL);
            xr = var(j).noturn.x_conv_lvl(rL);     
            lx = linspace(mean(xr)-3*std(xr), mean(xr)+3*std(xr), opstruct.numLxBins);
            Rates.rfLowVar_lvl(j)  = turnRates (xt, xr, deltaT, lx, opstruct);
        end
        
%        xt = var(j).turn.x_conv_ber120(tH);
%        xr = var(j).noturn.x_conv_ber120(rH);
%        lx = linspace(-3*std(xr), 3*std(xr), opstruct.numLxBins);   
%        Rates.rfHighVar_ber120(j) = turnRates (xt, xr, deltaT, lx, opstruct);
        
        xt = var(j).turn.x_conv(tL);
        xr = var(j).noturn.x_conv(rL);     
        lx = linspace(-3*std(xr), 3*std(xr), opstruct.numLxBins);
        Rates.rfLowVar(j)  = turnRates (xt, xr, deltaT, lx, opstruct);
        
        xt = var(j).turn.x_conv_LK(tL);
        xr = var(j).noturn.x_conv_LK(rL);     
        lx = linspace(-3*std(xr), 3*std(xr), opstruct.numLxBins);
        Rates.rfLowVar_LK(j)  = turnRates (xt, xr, deltaT, lx, opstruct);
                
%        xt = var(j).turn.x_conv_ber120(tL);
%        xr = var(j).noturn.x_conv_ber120(rL);     
%        lx = linspace(-3*std(xr), 3*std(xr), opstruct.numLxBins);
%        Rates.rfLowVar_ber120(j)  = turnRates (xt, xr, deltaT, lx, opstruct);
        
        if(opstruct.adaptationTime > 0)
            xt = var(j).turn.x_conv(tU);
            xr = var(j).noturn.x_conv(rU);
            lx = linspace(-3*std(xr), 3*std(xr), opstruct.numLxBins);
            Rates.rfUp(j) = turnRates (xt, xr, deltaT, lx, opstruct);
            
            xt = var(j).turn.x_conv(tD);
            xr = var(j).noturn.x_conv(rD);
            lx = linspace(-3*std(xr), 3*std(xr), opstruct.numLxBins);
            Rates.rfDown(j) = turnRates (xt, xr, deltaT, lx, opstruct);
            
        end
    end
    
    if(size(var, 2) == 2)    
        xt2D(:,1) = var(1).turn.x_conv;
        xt2D(:,2) = var(2).turn.x_conv;
        xr2D(:,1) = var(1).noturn.x_conv;
        xr2D(:,2) = var(2).noturn.x_conv;
        Rates.Rates2D.rfHighVar = turnRates2D (xt2D, xr2D, deltaT, tH, rH, opstruct.numLxBins);
        Rates.Rates2D.rfLowVar = turnRates2D (xt2D, xr2D, deltaT, tL, rL, opstruct.numLxBins);
        if(opstruct.adaptationTime > 0)
            Rates.Rates2D.rfUp = turnRates2D (xt2D, xr2D, deltaT, tU, rU, opstruct.numLxBins);
            Rates.Rates2D.rfDown = turnRates2D (xt2D, xr2D, deltaT, tD, rD, opstruct.numLxBins);
        end
    end
    
    btdstruct.Rates = Rates;
    
else
    
    for j = 1:size(var, 2)
        deltaT = median(diff(var(j).fullensemble.eti));
        
        xt = var(j).turn.x_conv;
        expnt = var(j).turn.expnum;
        xr = var(j).noturn.x_conv;
        expnr = var(j).noturn.expnum;
                
        lx = linspace(-3*std(xr), 3*std(xr), opstruct.numLxBins);
        
        Rates.rf(j) = turnRates (xt, xr, deltaT, lx, opstruct);
        
        for k=1:max(expnt)
           
            xt_sep = xt(expnt == k);
            xr_sep = xr(expnr == k);
            lx = linspace(-3*std(xr_sep), 3*std(xr_sep), opstruct.numLxBins);
            RatesSepExp{k}.rf(j) = turnRates (xt_sep, xr_sep, deltaT, lx, opstruct);
        end
        
    end
    btdstruct.RatesSepExp = RatesSepExp;
    btdstruct.Rates = Rates;
    
end

end

function [btdstruct, rateVtime] = measureRateVtime(var, deltaT, dt)

if( isfield(var.turn, 'ton') )
    tton = var.turn.ton;
    aton = var.noturn.ton;
    tx_ton = 0:dt:var.period;
    rateVtime.rate_ton = 60* histc(tton, tx_ton)./histc(aton, tx_ton) / deltaT;
    rateVtime.tx_ton = tx_ton;
end

teti = var.turn.eti;
aeti = var.noturn.eti;
tx = 0:dt:max(aeti);
rateVtime.rate = 60* histc(teti, tx)./histc(aeti, tx) / deltaT;

rateVtime.tx = tx;
btdstruct.rateVtime = rateVtime;

end

function [meanTTA, meanBS, stdBS] = bootstrapTTAfromBTDM(BTDM,expnum,tnum,taxis,numBootSamp)


datamat = BTDM; 
datamat = datamat';  %(assumes BTDM is in the form (turns x taxis) )

expmat = repmat(expnum,[1 length(taxis)]);  %matrix of experiment numbers
expmat = expmat';

tnummat = repmat(tnum,[1 length(taxis)]);  % matrix of track numbers
tnummat = tnummat';

%numtracks = cell(length(numBootSamp*max(expnum)),1);
drawtracks = cell(length(numBootSamp*max(expnum)),1);
for i=1:numBootSamp*max(expnum)   
    drawexp(i) = round( (max(expnum)-1)*rand) +1;  %array of random experiments to sample
    numtracks(i) = max( tnum(expnum == drawexp(i) ));  %number of tracks in each experiment in drawexp
    drawtracks{i} = round((numtracks(i)-1)*rand(1,numtracks(i)))+1;  %cell array of random tracks to draw from each experiment in drawexp
end

    
for k = 1:numBootSamp
    clear sampledat;
    sampledat = datamat(expmat == drawexp(1) & tnummat == drawtracks{1}(1));
    for l = 2:length(drawtracks{1}) 
        sampledat = cat(1,sampledat, datamat(expmat == drawexp(1) & tnummat == drawtracks{1}(l)));
    end
    for i=2:max(expnum)  
        for j=1:length(drawtracks{(k-1)*max(expnum)+i})
            sampledat = cat(1, sampledat, datamat(expmat == drawexp((k-1)*max(expnum)+i) & tnummat == drawtracks{(k-1)*max(expnum)+i}(j)));
        end
    end
    sampledat(~isfinite(sampledat)) = 0;
    sdatmat = vec2mat(sampledat,length(taxis));
    meanTTA(k,:) = mean(sdatmat,1);
%    samdatmat{k} = sdatmat;
end

meanBS = mean(meanTTA, 1);
stdBS = std(meanTTA , [], 1);

end