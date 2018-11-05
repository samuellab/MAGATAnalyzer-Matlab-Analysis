function [outputStruct] = fitRateFunWithTemporalScaling (ratefun, dlogratefun, ddlogratefun, turnT, turnX, runT, runX,  tx, deltaT, Q_alpha, params_0, alpha_0, w_0)
%function [alpha, valpha, staticRateParams, outputStruct] = fitRateFunWithTemporalScaling (ratefun, dlogratefun,ddlogratefun, turnT, turnX, runT, runX,  tx, deltaT, Q_alpha, params_0, alpha_0, w_0)
%
%fitRateFunWithTemporalScaling (ratefun, dlogratefun, ddlogratefun, turnT, turnX, runT, runX, tx, deltaT, Q_alpha, params_0, alpha_0, w_0)
%fits rate = ratefun(staticRateParams, xdata*alpha(t))
%logratefun - log (lambda(params, xdata))
%dlogratefun - d/dxdata (log (lambda(params, xdata))
%ddlogratefun - d^2/dxdata^2 (log (lambda(params, xdata))
%turnT, turnX, time and value at turn start
%runT, runX, time and value at runs (but not turns)
%tx time axis for analysis; data is assumed to be periodic with T =
%max(tx)-min(tx)
%deltaT - sampling interval for data (not related directly to tx) - usually 1/20
%Q_alpha - diffusion constant for scaling constant
%params0 - initial guess for rate function params
%slpha0 - initial guess for alpha(t)
%
%implements stochastic point process filter from eden et al 2004
if (nargin == 0)
    pd.ratefun = @(x, xdata) exp(polyval(x, xdata));
    pd.dlogratefun =@(x,xdata) polyval(polyder(x), xdata);
    pd.ddlogratefun =@(x,xdata) polyval(polyder(polyder(x)), xdata);
    pd.turnT = [];
    pd.turnX = [];
    pd.runT = [];
    pd.runX = [];
    pd.tx = [];
    pd.deltaT = 1/20;
    pd.Q_alpha = 1e-3;
    pd.params_0 = [0 1];
    pd.alpha_0 = ones(size(pd.tx));
    pd.w_0 = 1E9;
    pd.maxreps = 10;
    pd.pad = true;
    outputStruct = pd;
    return;
end

pad = true;
maxreps = 10;

if (nargin == 1 && isstruct(ratefun))
    pd = ratefun;
    fn = fieldnames(pd);
    for j = 1:length(fn);
        eval([fn{j} ' = pd.' fn{j} ';']);
    end
end

params = params_0;
alpha = alpha_0;
%llold = -Inf;
ll = NaN([1 maxreps]);
T = max(tx) - min(tx);
tinds = turnT >= min(tx) & turnT <= max(tx);
if (pad)
    turnXe = [turnX(tinds); turnX(tinds); turnX(tinds)];
    turnTe = [turnT(tinds)-T; turnT(tinds); turnT(tinds) + T];
    rinds = runT >= min(tx) & runT <= max(tx);
    runXe = [runX(rinds); runX(rinds); runX(rinds)];
    runTe = [runT(rinds)-T; runT(rinds);runT(rinds)+T];
    txe = unique([tx-T;tx;tx+T]);
else
    turnXe = turnX(tinds);
    turnTe = turnT(tinds);
    rinds = runT >= min(tx) & runT <= max(tx);
    runXe = runX(rinds);
    runTe = runT(rinds);
    txe = unique(tx);
end
    
[turnTe, I] = sort(turnTe, 'ascend');
turnXe = turnXe(I);
[runTe, I] = sort(runTe, 'ascend');
runXe = runXe(I);

[turnT, I] = sort(turnT, 'ascend');
turnX = turnX(I);
[runT, I] = sort(runT, 'ascend');
runX = runX(I);


for mm = 1:maxreps
    alpha = alpha /mean(alpha);
    [params, ll(mm)] = fitStaticRate (ratefun, params, turnT, turnX, runT, runX, tx, deltaT, alpha);
    fitstruct = fitScalingFactor(ratefun, dlogratefun, ddlogratefun, turnTe, turnXe, runTe, runXe, deltaT, txe, Q_alpha, params, alpha(1), w_0);
    alpha = interp1(txe, fitstruct.theta, tx, 'nearest');
%     alpha = alpha /mean(alpha);
    w_0 = interp1(txe, fitstruct.w, tx(1), 'nearest');
    plot (ll); pause(0.01);
end
%'tx','alpha', 'valpha', 'alpha_s', 'valpha_s',
fitstruct.tx = tx;
fitstruct.ll = ll;
fitstruct.alpha = interp1(txe, fitstruct.theta, tx, 'nearest'); 
fitstruct.alpha_s = interp1(txe, fitstruct.theta_s, tx, 'nearest');
fitstruct.valpha = interp1(txe, fitstruct.w, tx, 'nearest');
fitstruct.valpha_s = interp1(txe, fitstruct.w_s, tx, 'nearest');

outputStruct = fitstruct;
outputStruct.staticParams = params;






function [params, ll] = fitStaticRate (ratefun, params_0, turnT, turnX, runT, runX, tx, deltaT, alpha)
%function params = fitStaticRate (logratefun, params_0, turnT, turnX, runT, runX, tx, alpha)
%maximizes log(P(data|params)) = sum_turn logratefun(params, turnX*alpha) - sum_run exp(logratefun(params,runX*alpha)))*deltaT 
tval = turnX(turnT >= min(tx) & turnT <= max(tx)).*interp1(tx, alpha, turnT(turnT >= min(tx) & turnT <= max(tx)), 'linear');
rval = runX(runT >= min(tx) & runT <= max(tx)).*interp1(tx, alpha, runT(runT >= min(tx) & runT <= max(tx)), 'linear');
nlogP = @(p) -sum(log(ratefun(p, tval)) ) + sum(ratefun(p,rval)*deltaT);

op = optimoptions('fminunc');
op.Display = 'off';
op.Algorithm = 'quasi-newton';

[params, ll, exitflag] = fminunc(nlogP, params_0, op);
% nt = length(tval)
% ntpred = sum(ratefun(params,rval)*deltaT)
% innov = zeros(size(tx));
% nt = innov;
% ntpred = innov;
% % [~,~,tbin] = histcounts(turnT, tx);
% % [~,~,rbin] = histcounts(runT, tx);
% 
% 
% for j = 1:(length(tx) - 1)
%     ti = findRangeInSortedData(turnT, tx(j), tx(j+1));
%     ri = findRangeInSortedData(runT, tx(j), tx(j+1));
% %     td = turnX(ti);
%     ad = runX(ri);
% %     xt = td*alpha(j);
%     xa = ad*alpha(j);
%     nt(j) =  length(ti);
%     ntpred(j) = deltaT*sum(ratefun(params, xa));
% end
% innov = nt - ntpred;
% mean(innov)
% plot (innov);


ll = -ll + length(tval)*log(deltaT);
if (exitflag <= 0)
    warning ('static rate fit may not have converged!');
end


function fitstruct = fitScalingFactor (ratefun, dlogratefun,ddlogratefun, turnT, turnX, runT, runX, deltaT, txe, Q_alpha, params, alpha_0, w_0)
existsAndDefault('alpha_0', 1);


% allT = [turnT runT];
% allX = [turnX runX];
% isturn = [true(size(turnT)) false(size(runT))];
theta = ones(size(txe));
w = theta;
wkkm1 = w;
innov = zeros(size(theta));


theta(1) = alpha_0;
w(1) = w_0;
wkkm1(1) = w_0;
% [~,~,tbin] = histcounts(turnT, txe);
% [~,~,rbin] = histcounts(runT, txe);

for j = 1:(length(txe) - 1)
%    t = theta(j);
    wkkm1(j+1) = w(j) + Q_alpha*(txe(j+1)-txe(j));  
    ti = findRangeInSortedData(turnT, txe(j), txe(j+1));
    ri = findRangeInSortedData(runT, txe(j), txe(j+1));
%     ti = tbin == j;
%     ri = rbin == j;
    td = turnX(ti);
%     ad = [td;runX(ri)];
    ad = runX(ri);
    xt = td*theta(j);
    xa = ad*theta(j);
    w(j+1) = 1 / (1/wkkm1(j+1) ...
        + deltaT*sum(ratefun(params, xa).*((dlogratefun(params, xa)).^2 + ddlogratefun(params, xa)).*ad.^2) ...
        - sum(ddlogratefun(params, xt).*td.^2) );
    theta(j+1) = theta(j) + w(j+1)*(sum(dlogratefun(params, xt).*td) - deltaT*sum(ratefun(params, xa).*dlogratefun(params, xa).*ad));
    innov(j) = length(ti) - deltaT*sum(ratefun(params, xa));
    
end

[theta_s,w_s] = recursiveSmoothing(theta, w, wkkm1);

fn = { 'txe', 'theta', 'w', 'innov', 'wkkm1', 'theta_s', 'w_s'};
for j = 1:length(fn)
    fitstruct.(fn{j}) = eval(fn{j}); 
end




function [ts, ws] = recursiveSmoothing (t, w, wkkm1)
%function [ts, ws] = recursiveSmoothing (t, w, wkkm1)
%
%implements smoothing step from Koyama, S., Eden, U.T., Brown, E.N., and Kass, R.E. (2009). Bayesian decoding of neural spike trains. Ann Inst Stat Math 62, 37.
%eqn 17-19
%assuming f = 1
ts = t;
ws = w;

for j = (length(w)-1):-1:1
    h = w(j)/wkkm1(j+1);
    ts(j) = t(j) + h*(ts(j+1) - t(j));
    ws(j) = ws(j) + h^2*(ws(j+1)-wkkm1(j+1));
end


