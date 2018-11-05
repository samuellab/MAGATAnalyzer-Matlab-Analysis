function [outputStruct, scaledValues] = fitRateFunWithTemporalScalingND (pd)
%function [outputStruct] = fitRateFunWithTemporalScalingND (pd)
%
% pd.ratefun = @(x, xdata) exp(x(1)*xdata(:,1) + x(2)*xdata(:,2) + x(3)); %[Nx1 rate function of NxD xdata given K params]
% pd.gradlogratefun =@(x,xdata) [x(1)*ones(size(xdata(:,1))) x(2)*ones(size(xdata(:,2)))] ; %NxD gradient of log of rate
% pd.hesslogratefun = @(x,xdata) zeros([size(xdata,1) size(xdata,2) size(xdata,2)]); %NxDxD hessian of log of rate;
% pd.temporalratemod = @(x,tdata) x(1)*exp(-tdata/x(2)) + 1;
% pd.turnT = [];
% pd.turnEti = [];
% pd.turnX = [];
% pd.runT = [];
% pd.runEti = [];
% pd.runX = [];
% pd.tx = [];
% pd.etx = [];
% pd.deltaT = 1/20;
% pd.Q_alpha = 1e-3;
% pd.params_0 = [0 0 1];
% pd.tparams_0 = [0 1000];
% pd.alpha_0 = ones(size(pd.tx));
% pd.w_0 = 1E9*ones(size(pd.turnX,2)*[1 1]);
% pd.maxreps = 10;
% pd.pad = false;
% outputStruct = pd;
% pd.separateExperiments = false;
%
%implements stochastic point process filter from eden et al 2004
if (nargin == 0)
    pd.ratefun = @(x, xdata) exp(x(1)*xdata(:,1) + x(2)*xdata(:,2) + x(3)); %[Nx1 rate function of NxD xdata given K params]
    pd.gradlogratefun =@(x,xdata) [x(1)*ones(size(xdata(:,1))) x(2)*ones(size(xdata(:,2)))] ; %NxD gradient of log of rate
    pd.hesslogratefun = @(x,xdata) zeros([size(xdata,1) size(xdata,2) size(xdata,2)]); %NxDxD hessian of log of rate;
    pd.temporalratemod = @(x,tdata) x(1)*exp(-tdata/x(2)) + 1;
    pd.turnT = [];
    pd.turnEti = [];
    pd.turnX = [];
    pd.runT = [];
    pd.runEti = [];
    pd.runX = [];
    pd.tx = [];
    pd.etx = [];
    pd.deltaT = 1/20;
    pd.Q_alpha = 1e-3;
    pd.params_0 = [0 0 1];
    pd.tparams_0 = [0 1000];
    pd.alpha_0 = ones(size(pd.tx));
    pd.w_0 = 1E9*ones(size(pd.turnX,2)*[1 1]);
    pd.maxreps = 10;
    pd.pad = false;
    outputStruct = pd;
    pd.separateExperiments = false;
    return;
end
% 
% pad = true;
% maxreps = 10;

% if (nargin == 1 && isstruct(ratefun))
%     pd = ratefun;
separateExperiments = false;
fn = fieldnames(pd);
for j = 1:length(fn)
    eval([fn{j} ' = pd.' fn{j} ';']);
end
% end

if isempty(runEti) %#ok<*NODEF>
    runEti = runT;
end
if isempty(turnEti)
    turnEti = turnT;
end

tx = tx(:);

%params = params_0;
alpha = alpha_0;
%llold = -Inf;
ll = NaN([1 maxreps]);
T = max(tx) - min(tx);
tinds = turnT >= min(tx) & turnT < max(tx);
rinds = runT >= min(tx) & runT < max(tx);

turnT = turnT(tinds);
runT = runT(rinds);
turnEti = turnEti(tinds);
runEti = runEti(rinds);
turnX = turnX(tinds,:);
runX = runX(rinds,:);
turnExpnum = turnExpnum(tinds);
runExpnum = runExpnum(rinds);

[turnT, I] = sort(turnT, 'ascend');
turnEti = turnEti(I);
turnX = turnX(I, :);
turnExpnum = turnExpnum(I);
[runT, I] = sort(runT, 'ascend');
runEti = runEti(I);
runX = runX(I, :);
runExpnum = runExpnum(I);



if (pd.pad)
    turnXe = [turnX; turnX; turnX];
    turnEtie = [turnEti;turnEti;turnEti];
    turnTe = [turnT-T; turnT; turnT + T];  
    runXe = [runX; runX; runX];
    runTe = [runT-T; runT;runT+T];
    runEtie = [runEti;runEti;runEti];
    turnExpnume = [turnExpnum; turnExpnum; turnExpnum];
    runExpnume = [runExpnum; runExpnum; runExpnum];
    txe = unique([tx-T;tx;tx+T]);
else
    turnXe = turnX;
    turnTe = turnT;
    turnEtie = turnEti;
    runXe = runX;
    runTe = runT;
    runEtie = runEti;
    turnExpnume = turnExpnum; 
    runExpnume = runExpnum; 
    txe = unique(tx);
end

fn = {'turnT', 'turnX', 'turnEti', 'runT', 'runX', 'runEti', 'tx', 'turnExpnum', 'runExpnum'};
for j = 1:length(fn)
    data.(fn{j}) = eval(fn{j});
    data.([fn{j} 'e']) = eval([fn{j} 'e']);
end
data.deltaT = deltaT;

    
funs.ratefun = ratefun;
funs.temporalratemod = temporalratemod;
funs.params = params_0;
funs.tparams = tparams_0;
funs.gradlogratefun = gradlogratefun; %NxD gradient of log of rate
funs.hesslogratefun = hesslogratefun; %DxDxN hessian of log of rate;

if(size(alpha, 2)>2)
    alpha = alpha';
end
theta_0 = alpha(1,:);
nr_t = histcounts(pd.runT, binEdgesFromCenters(pd.tx));
nt_t = histcounts(pd.turnT, binEdgesFromCenters(pd.tx));
wt = (nr_t + nt_t)/sum(nr_t + nt_t);
wt = wt(:);
norm_factor = sum(alpha.*repmat(wt, [1 size(alpha,2)]),1);    
alpha = alpha ./ repmat(norm_factor, [size(alpha,1) 1]);

[data.tval, data.rval] = stretchedValues (data.tx, alpha, data.turnT, data.turnX, data.runT, data.runX);


nr_et = histcounts(pd.runEti, binEdgesFromCenters(pd.etx));
nt_et = histcounts(pd.turnEti, binEdgesFromCenters(pd.etx));
tr = 1/data.deltaT .* nt_et./nr_et;
if (~isempty(tparams_0) && ~isempty(temporalratemod))
    tratefun = @(p, td) p(1) * temporalratemod(p(2:end),td);
    tp = lsqcurvefit(tratefun, [mean(tr) tparams_0], pd.etx, tr);
    funs.tparams = tp(2:end);
end
tic
funs = fitStaticRate (funs, data);
toc
if (separateExperiments)
    erange = [min(data.turnExpnum) max(data.turnExpnum)];
else
    erange = 1;
end

for mm = 1:maxreps
        for j = min(erange):max(erange)
            if (mm > 1)
                w_0 = fitstruct(j).w_s(:,:,1);
                theta_0 = fitstruct(j).alpha_s(1,:);
            end
            if (separateExperiments)
                fitstruct(j) = fitScalingFactor(funs, data, Q_alpha, theta_0, w_0,j);  %#ok<AGROW>
            else
                fitstruct = fitScalingFactor(funs, data, Q_alpha, theta_0, w_0); 
            end
        end
     toc
     data = addStretchedValues(data, fitstruct);
    postAlphaFit = measuredAndPredictedTurns (funs, data, tx); 
%     
    fitstruct = normAlpha(fitstruct, data);
    data = addStretchedValues(data, fitstruct);
    [funs, ll(mm), nll(mm)] = fitStaticRate (funs, data);
    
     postRateFit = measuredAndPredictedTurns (funs, data, tx); 
     toc
%     
   
end
%'tx','alpha', 'valpha', 'alpha_s', 'valpha_s',
% fitstruct.tx = tx;
% fitstruct.ll = ll;

outputStruct.fitstruct = fitstruct;
outputStruct.staticParams = funs.params;
outputStruct.temporalParams = funs.tparams;
outputStruct.etx = etx;
outputStruct.turnRateVsTime = tr;
outputStruct.nll = nll;
outputStruct.ll = ll;
outputStruct.postAlphaFit = postAlphaFit;
outputStruct.postRateFit = postRateFit;

if (length(fitstruct) == 1)
    fn = fieldnames(fitstruct);
    for j = 1:length(fn)
        outputStruct.(fn{j}) = fitstruct.(fn{j});
    end
else
   tx = fitstruct(1).tx;   
   a = zeros([size(fitstruct(1).alpha) length(fitstruct) ]);
   a_s = a;
   v = zeros([size(fitstruct(1).valpha) length(fitstruct) ]);
   v_s = v;
   for j = 1:length(fitstruct)
       a(:,:,j) = fitstruct(j).alpha;
       a_s(:,:,j) = fitstruct(j).alpha_s;
       v(:,:,:,j) = fitstruct(j).valpha;
       v_s(:,:,:,j) = fitstruct(j).valpha_s;
   end
   alpha = fitstruct(1).alpha;
   alpha_s = alpha;
   valpha = fitstruct(1).valpha;
   valpha_s = valpha;
   for j = 1:length(tx)
       u = squeeze (a(j,:,:));
       if(size(u,2) == 1) %squeezed too far
           u = u';
       end
       wi = zeros([size(v,1) size(v,2), size(v,4)]);     
       for k = 1:size(v,4)
           wi(:,:,k) = inv(v(:,:,j,k));
           u(:,k) = wi(:,:,k)*u(:,k);
       end
       valpha(:,:,j) = inv(sum(wi,3));
       alpha(j,:) = (valpha(:,:,j)*sum(u,2))';
       
       u = squeeze (a_s(j,:,:));
       if(size(u,2) == 1) %squeezed too far
           u = u';
       end
       wi = zeros([size(v,1) size(v,2), size(v,4)]);     
       for k = 1:size(v,4)
           wi(:,:,k) = inv(v_s(:,:,j,k));
           u(:,k) = wi(:,:,k)*u(:,k);
       end
       valpha_s(:,:,j) = inv(sum(wi,3));
       alpha_s(j,:) = (valpha_s(:,:,j)*sum(u,2))';
   end
   
   outputStruct.tx = tx;
   outputStruct.innov = sum([fitstruct.innov], 2);     
   outputStruct.alpha = alpha;
   outputStruct.alpha_s = alpha_s;
   outputStruct.valpha = alpha;
   outputStruct.valpha_s = alpha_s;
   outputStruct.pd = pd;
   
   
end
scaledValues = data;




function [funs, ll, nll] = fitStaticRate (funs, data)
%function params = fitStaticRate (logratefun, params_0, turnT, turnX, runT, runX, tx, alpha)
%maximizes log(P(data|params)) = sum_turn logratefun(params, turnX*alpha) - sum_run exp(logratefun(params,runX*alpha)))*deltaT 

tval = data.tval;
rval = data.rval;
teti = data.turnEti;
reti = data.runEti;

op = optimoptions('fminunc');
op.Display = 'off';
op.Algorithm = 'quasi-newton';



if (isempty(funs.temporalratemod) || isempty(funs.tparams))
    funs.tparams = [];
    funs.temporalratemod = @(p,td) ones(size(td));
    ratefun = @(p,xd, td) funs.ratefun(p, xd);
    p0 = funs.params;
    
    nlogP = @(p) -sum(log(ratefun(p, tval, teti)) ) + sum(ratefun(p,rval, reti)*data.deltaT);
    [funs.params, ll, exitflag] = fminunc(nlogP, p0, op);
else
    tparams = funs.tparams;
    params = funs.params;
    ratefun = @(p,tp,xd, td) funs.ratefun(p, xd) .* funs.temporalratemod(tp, td);
    nlogP = @(p)  -sum(log(ratefun(p, tparams, tval, teti)) ) + sum(ratefun(p,tparams, rval, reti)*data.deltaT);
    [params, ~, exitflag] = fminunc(nlogP, params, op);
    if (exitflag <= 0)
        warning ('static rate fit may not have converged!');
    end
    nlogP = @(tp)  -sum(log(ratefun(params, tp, tval, teti)) ) + sum(ratefun(params,tp, rval, reti)*data.deltaT);
    [tparams, ~, exitflag] = fminunc(nlogP, tparams, op);
    if (exitflag <= 0)
        warning ('static rate fit may not have converged!');
    end
    nlogP = @(p)  -sum(log(ratefun(p(1:length(params)), p((length(params)+1):end), tval, teti)) ) + sum(ratefun(p(1:length(params)), p((length(params)+1):end), rval, reti)*data.deltaT);
    p0 = [params tparams];
    [p, ll, exitflag] = fminunc(nlogP, p0, op);    
    funs.params = p(1:length(funs.params));
    funs.tparams = p((length(funs.params) + 1):end);
    
end


nll = ll;
ll = -ll + length(tval)*log(data.deltaT);
if (exitflag <= 0)
    warning ('static rate fit may not have converged!');
end


function fitstruct = fitScalingFactor(funs, data, Q_alpha, alpha_0, w_0, expnum) % (ratefun, dlogratefun,ddlogratefun, turnT, turnX, runT, runX, deltaT, txe, Q_alpha, params, alpha_0, w_0)
existsAndDefault('alpha_0', 1);


txe = data.txe;
theta = repmat(alpha_0, [length(txe), 1]);
dtheta = theta;
%theta = ones(size(txe));
w = repmat(w_0, [1 1 length(txe)]);
wkkm1 = w;
innov = zeros(length(txe), 1);

if (numel(Q_alpha) < numel(w_0))
    Q_alpha = Q_alpha(1) * eye(size(w_0,1));
end

if (existsAndDefault('expnum', []))
    tinds = data.turnExpnume == expnum;
    rinds = data.runExpnume == expnum;
else
    tinds = true(size(data.turnTe));
    rinds = true(size(data.runTe));
end
turnTe = data.turnTe(tinds);
turnEtie = data.turnTe(tinds);
turnXe = data.turnXe(tinds, :);


runTe = data.runTe(rinds);
runEtie = data.runTe(rinds);
runXe = data.runXe(rinds, :);


for j = 1:(length(txe) - 1)
%    t = theta(j);
    wkkm1(:,:,j+1) = w(:,:,j) + Q_alpha*(txe(j+1)-txe(j));  
    ti = findRangeInSortedData(turnTe, txe(j), txe(j+1));
    ri = findRangeInSortedData(runTe, txe(j), txe(j+1));
    if (isempty(ri))
        w(:,:,j+1) = wkkm1(:,:,j+1);
        theta(j+1,:) = theta(j,:);
        continue;
    end
    
    td = turnXe(ti, :);
    rd = runXe(ri, :);
    %ad = runXe(ri);
    ad = [td;rd];
    xt = td.*repmat(theta(j, :), [size(td,1) 1]);
    xr = rd.*repmat(theta(j, :), [size(rd,1) 1]);
    xa = ad.*repmat(theta(j, :), [size(ad,1) 1]);
    
    teti = turnEtie(ti);
    reti = runEtie(ri);
    aeti = [teti;reti];
    
    %tlambda = funs.ratefun(funs.params, xt).*funs.temporalratemod(funs.tparams, teti);
    rlambda_dt = data.deltaT*funs.ratefun(funs.params, xr).*funs.temporalratemod(funs.tparams, reti);
    alambda_dt = data.deltaT*funs.ratefun(funs.params, xa).*funs.temporalratemod(funs.tparams, aeti);
    innov(j) = length(ti) - sum(rlambda_dt);
    
    
    ahess = funs.hesslogratefun(funs.params, xa);
    if (isempty(ahess))
        ahesstheta = ahess;
    else
        adm = repmat(ad, [1 1 size(ad,2)]);
        ahesstheta = ahess.*adm.*permute(adm, [1 3 2]); %Hij * xi * xj 
    end
    thess = funs.hesslogratefun(funs.params, xt);
    if (isempty(thess))
        thesstheta = thess;
    else
        tdm = repmat(td, [1 1 size(td,2)]);
        thesstheta = thess.*tdm.*permute(tdm, [1 3 2]);
    end
    
    agrad = funs.gradlogratefun(funs.params, xa);
    if (isempty(agrad))
        agradthetasq = ahesstheta;
    else
        agradtheta = agrad.*ad;
        agradthetasq = repmat(agradtheta, [1 1 size(agradtheta,2)]);
        agradthetasq = agradthetasq.*permute(agradthetasq, [1 3 2]);
    end
    
%     dlambdaA_lambdaA = zeros(size(agradtheta,2)*[1 1]);
%     for m = 1:size(agradtheta,2)
%         for n = m:size(agradtheta,2)
%             dlambdaA_lambdaA(m,n) = sum(agradtheta(:,m).*agradtheta(:,n).*alambda_dt);
%             if (m < n)
%                 dlambdaA_lambdaA(n,m) = dlambdaA_lambdaA(m,n);
%             end
%             
%         end
%     end
%     dlambdaA_lambdaA_rev = squeeze(sum((agradthetasq).*repmat(alambda_dt, [1 size(ahesstheta, 2) size(ahesstheta,2)]),1));
%     
    tgrad = funs.gradlogratefun(funs.params, xt);
    tgradtheta = tgrad.*td;
    
    rgrad = funs.gradlogratefun(funs.params, xr);
    rgradtheta = rgrad.*rd;
    
    
    dwinv = squeeze(sum((agradthetasq + ahesstheta).*repmat(alambda_dt, [1 size(ahesstheta, 2) size(ahesstheta,2)]),1) - ...
        sum(thesstheta,1));
    
    w(:,:,j+1) = inv(inv(wkkm1(:,:,j+1)) + dwinv);
    dtheta(j,:) = (w(:,:,j+1)*(sum(tgradtheta,1)-sum(rgradtheta .* repmat(rlambda_dt, [1 size(rgradtheta, 2)]),1))')';
    theta(j+1,:) = theta(j,:) + dtheta(j,:); 
    
    
end

[theta_s,w_s] = recursiveSmoothing(theta, w, wkkm1); %#ok<ASGLU>

fn = { 'txe', 'theta', 'w', 'innov', 'wkkm1', 'theta_s', 'w_s', 'dtheta'};
for j = 1:length(fn)
    fitstruct.(fn{j}) = eval(fn{j}); 
end
fitstruct.tx = data.tx;
fitstruct.alpha = interp1(txe, fitstruct.theta, data.tx, 'nearest'); 
fitstruct.alpha_s = interp1(txe, fitstruct.theta_s, data.tx, 'nearest');
for j = 1:size(fitstruct.w, 1)
    for k = 1:size(fitstruct.w, 2)        
        fitstruct.valpha(j,k,:) = interp1(txe, squeeze(fitstruct.w(j,k,:)), data.tx, 'nearest');
        fitstruct.valpha_s(j,k,:) = interp1(txe, squeeze(fitstruct.w_s(j,k,:)), data.tx, 'nearest');
    end
end

function result = measuredAndPredictedTurns (funs, data, tx)

ntpred = zeros(size(tx));
nt = ntpred;
nr = ntpred;
for j = 1:(length(tx) - 1)
    ti = findRangeInSortedData(data.turnT, tx(j), tx(j+1));
    ri = findRangeInSortedData(data.runT, tx(j), tx(j+1));
    if (isempty(ri))
        continue;
    end    
    xa = data.rval(ri, :);
    %xa = ad.*repmat(alpha(j, :), [size(ad,1) min(size(ad,1),1)]);
    aeti = data.runEtie(ri);
    nt(j) = length(ti);
    nr(j) = length(ri);    
    ntpred(j) = sum(data.deltaT*funs.ratefun(funs.params, xa).*funs.temporalratemod(funs.tparams, aeti));
    
end

result.ntpred = ntpred;
result.nt = nt;
result.nr = nr;
result.innov = nt - ntpred;
result.residual = cumsum(result.innov);



function [ts, ws] = recursiveSmoothing (t, w, wkkm1)
%function [ts, ws] = recursiveSmoothing (t, w, wkkm1)
%
%implements smoothing step from Koyama, S., Eden, U.T., Brown, E.N., and Kass, R.E. (2009). Bayesian decoding of neural spike trains. Ann Inst Stat Math 62, 37.
%eqn 17-19
%assuming F = 1
ts = t;
ws = w;

for j = (length(w)-1):-1:1
    h = w(:,:,j)/wkkm1(:,:,j+1);
    ts(j,:) = t(j,:) + (h*(ts(j+1,:) - t(j,:))')';
    ws(:,:,j) = ws(:,:,j) + h*(ws(:,:,j+1)-wkkm1(:,:,j+1))*h';
end

function fitstruct = normAlpha (fitstruct, data, expnum)
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

nr_t = histcounts(data.runT(rinds), binEdgesFromCenters(data.tx));
nt_t = histcounts(data.turnT(tinds), binEdgesFromCenters(data.tx));
wt = (nr_t + nt_t)/sum(nr_t + nt_t);
wt = wt(:);
alpha = fitstruct.alpha;
norm_factor = sum(alpha.*repmat(wt, [1 size(alpha,2)]),1);
fitstruct.alpha = alpha ./ repmat(norm_factor, [size(alpha,1) 1]);

alpha_s = fitstruct.alpha_s;
norm_factor = sum(alpha_s.*repmat(wt, [1 size(alpha,2)]),1);
fitstruct.alpha_s = alpha_s ./ repmat(norm_factor, [size(alpha,1) 1]);


function data = addStretchedValues (data, fitstruct)
if (length(fitstruct) == 1)
    [data.tval, data.rval] = stretchedValues (fitstruct.tx, fitstruct.alpha, data.turnT, data.turnX, data.runT, data.runX);
else
    data.tval = data.turnX;
    data.rval = data.runX;
    for j = 1:length(fitstruct)
        ti = data.turnExpnum == j;
        ri = data.runExpnum == j;
        [tval, rval] = stretchedValues (fitstruct(j).tx, fitstruct(j).alpha, data.turnT(ti), data.turnX(ti,:), data.runT(ri), data.runX(ri,:));
        data.tval(ti,:) = tval;
        data.rval(ri,:) = rval;
    end
end

function [tval, rval] = stretchedValues (tx, alpha, turnT, turnX, runT, runX)
tval = turnX .* interp1(tx, alpha, turnT, 'linear');
rval = runX .* interp1(tx, alpha, runT, 'linear');
