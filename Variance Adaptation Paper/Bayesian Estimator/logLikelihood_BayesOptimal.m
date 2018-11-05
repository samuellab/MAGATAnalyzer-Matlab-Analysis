function SF = logLikelihood_BayesOptimal(btdstruct, dts, taus, bayes_sig, sigtoalpha, alphastruct)
%function SF = logLikelihood_BayesOptimal(btdstruct, dts, taus, bayes_sig, sigtoalpha, alphastruct)

deltaT = median(diff(btdstruct.var.fullensemble.eti));

SimData.sim = bayes_sig;
SimData.dts = dts;
SimData.taus = taus;
SimData.sigtoalpha = sigtoalpha;



tton = btdstruct.var.turn.ton';
rton = btdstruct.var.noturn.ton';
teti = btdstruct.var.turn.eti';
reti = btdstruct.var.noturn.eti';
xt = btdstruct.var.turn.x_conv';
xr = btdstruct.var.noturn.x_conv';

Data.alphastruct = alphastruct;

Data.tton = tton;
Data.rton = rton;
Data.xt = xt;
Data.xr = xr;

[fP, LL] = LL_fromFitData(Data, SimData, deltaT);
SF.dts = dts;
SF.taus = taus;
SF.fP = fP;
SF.LL = LL;

end

function [fP,LL] = LL_fromFitData(Data, SimData, deltaT)

op = optimoptions('fminunc');
op.Algorithm = 'quasi-newton';

sim = SimData.sim;
dts = SimData.dts;
taus = SimData.taus;
sigtoalpha = SimData.sigtoalpha;

alphastruct = Data.alphastruct;

gainFun = @(x, xdata) 1./(sqrt(xdata.^2 + x.^2))./mean(1./(sqrt(xdata.^2 + x.^2)));   
logR  = @(fitparams, xdata) polyval(fitparams, xdata);


for n=1:size(Data.tton, 1)
    tic
    
    teti = alphastruct.scaledValues.turnEti';
    reti = alphastruct.scaledValues.runEti';
    tton = alphastruct.scaledValues.turnTon';
    rton = alphastruct.scaledValues.runTon';
    xt = alphastruct.scaledValues.turnX';
    xr = alphastruct.scaledValues.runX';
    
    tx = alphastruct.tx;
    alpha = alphastruct.alpha;
    tval = xt.* interp1(tx, alpha, teti, 'spline');
    rval = xr.* interp1(tx, alpha, reti, 'spline');
    logPdata = @(params) -sum(logR(params, tval)) + deltaT*sum(exp(logR(params, rval)));
    fp0 = alphastruct.params;
    [fP.fp_data, logP_data] = fminunc(logPdata, fp0, op);
    
    tx_ton = alphastruct.tx_ton;
    alpha_ton = alphastruct.alpha_ton;
    tval = xt.* interp1(tx_ton, alpha_ton, tton, 'spline');
    rval = xr.* interp1(tx_ton, alpha_ton, rton, 'spline');
    logPdata = @(params) -sum(logR(params, tval)) + deltaT*sum(exp(logR(params, rval)));
%     [fP.fp_data_ton, logP_data_ton] = fminunc(logPdata, fp0, op);
    
    
    for i=1:length(dts)
        disp(['dt = ' num2str(dts(i))]);
        for k=1:length(taus)
            disp(['tau = ' num2str(taus(k))]);
            tic
            
            tx = sim{i}{k}.tx;
            inds = tx>=min(reti);
            alphaSim.tx = tx(inds);
            sigma = sim{i}{k}.sigma;
            alpha = @(StoA) gainFun(StoA, sigma(inds));
            
            %fit only rate params
            alphaSim.alpha = alpha(sigtoalpha);
            alphaSim = normAlpha(alphaSim, alphastruct.scaledValues, 'eti', [], [], []);
            tval = xt.* interp1(alphaSim.tx, alphaSim.alpha, teti, 'spline');
            rval = xr.* interp1(alphaSim.tx, alphaSim.alpha, reti, 'spline');
            logPdata = @(params) -sum(logR(params, tval)) + deltaT*sum(exp(logR(params, rval)));
            [fP.fp_fit(i,k,:), logP_fit(i,k)] = fminunc(logPdata, fP.fp_data, op);
            logP_nofit(i,k) = logPdata(fP.fp_data);
             
            toc
            %fit also alpha

%             findNormAlpha = @(s0) normAlpha(alphaSim, alphastruct.scaledValues, 'eti', [], sigma(inds), s0);
%             tval = @(s0) xt.* interp1(alphaSim.tx, findNormAlpha(s0), teti, 'spline');
%             rval = @(s0) xr.* interp1(alphaSim.tx, findNormAlpha(s0), reti, 'spline');
%             logPdata = @(params) -sum(logR(params(2:end), tval(params(1)))) + deltaT*sum(exp(logR(params(2:end), rval(params(1)))));
%             [fP.fp_fit(i,k,:), logP_fit_alpha(i,k)] = fminunc(logPdata, [1.5 fP.fp_data], op);
%             disp(['s0 = ' num2str(fP.fp_fit(i,k,1))]);
%             toc
%             logP_nofit(i,k) = logPdata(fP.fp_data);
            
%             alphaSim.tx_ton = sim{i}{k}.tx_ton;
%             sigma_ton = sim{i}{k}.sigma_ton;
%             alphaSim.alpha_ton = gainFun(sigtoalpha, sigma_ton);
%             alphaSim = normAlpha(alphaSim, alphastruct.scaledValues, 'ton', []);
%             tval = xt.* interp1(alphaSim.tx_ton, alphaSim.alpha_ton, tton, 'spline');
%             rval = xr.* interp1(alphaSim.tx_ton, alphaSim.alpha_ton, rton, 'spline');
%             logPdata = @(params) -sum(logR(params, tval)) + deltaT*sum(exp(logR(params, rval)));
%             [fP.fp_fit_ton(i,k,:), logP_fit_ton(i,k)] = fminunc(logPdata, fP.fp_data_ton, op);
%             logP_nofit_ton(i,k) = logPdata(fP.fp_data_ton);
            
         end
    end
end

LL.logP_data(n) = logP_data;
% LL.logP_data_ton(n) = logP_data_ton;

LL.logP_fitR(:,:,n) = logP_fit;
LL.logP_nofit(:,:,n) = logP_nofit;
% LL.logP_fitR_ton(:,:,n) = logP_fit_ton;
% LL.logP_nofit_ton(:,:,n) = logP_nofit_ton;

end





function fitstruct = normAlpha (fitstruct, data, timeType, expnum, sigma, s0)
if (length(fitstruct) > 1)
    for j = 1:length(fitstruct)
        fitstruct(j) = normAlpha (fitstruct(j), data, j);
    end
    return;
end
if (existsAndDefault('expnum', []))
    tinds = data.turnExpnum == expnum;
    rinds = data.runExpnum == expnum;
else
    tinds = true(size(data.turnT));
    rinds = true(size(data.runT));
end

if(~isfield(fitstruct, 'alpha'))
    gainFun = @(x) 1./(sqrt(sigma.^2 + x.^2))./mean(1./(sqrt(sigma.^2 + x.^2)));  
    fitstruct.alpha = gainFun(s0);
end

if(strcmpi(timeType, 'ton'))
    if(isfield(data, 'runTon'))
        runT = data.runTon;
        turnT = data.turnTon;
    else
        runT = data.runT;
        turnT = data.turnT;
    end
    nr_ton = histcounts(runT(rinds), binEdgesFromCenters(fitstruct.tx_ton));
    nt_ton = histcounts(turnT(tinds), binEdgesFromCenters(fitstruct.tx_ton));
    wton = (nr_ton + nt_ton)/sum(nr_ton + nt_ton);
    wton = wton(:);
    alpha_ton = fitstruct.alpha_ton';
    norm_factor = sum(alpha_ton.*repmat(wton, [1 size(alpha_ton,2)]),1);
    alpha_ton = alpha_ton ./ repmat(norm_factor, [size(alpha_ton,1) 1]);
    fitstruct.alpha_ton = alpha_ton';
else
    nr_t = histcounts(data.runT(rinds), binEdgesFromCenters(fitstruct.tx));
    nt_t = histcounts(data.turnT(tinds), binEdgesFromCenters(fitstruct.tx));
    wt = (nr_t + nt_t)/sum(nr_t + nt_t);
    wt = wt(:);
    alpha = fitstruct.alpha';
    norm_factor = sum(alpha.*repmat(wt, [1 size(alpha,2)]),1);
    alpha = alpha ./ repmat(norm_factor, [size(alpha,1) 1]);
    fitstruct.alpha = alpha';
end

if(~isempty(s0))
    fitstruct = fitstruct.alpha;
end


end

