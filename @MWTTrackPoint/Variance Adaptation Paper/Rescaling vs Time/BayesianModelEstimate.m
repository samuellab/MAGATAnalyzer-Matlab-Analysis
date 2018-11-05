
function out = BayesianModelEstimate(pd, btdstruct, var, opstruct, timeType, dt, Ddt, taxis, logdat, output, rescale_type, multiType)
% function out = BayesianModelEstimate(pd, btdstruct, var, opstruct, timeType, dt, Ddt, taxis, logdat, output, rescale_type)
% Calculates rescaling model parameter (alpha(t)) using a Bayesian estimator, fits the rate function to r(alpha(t),x), and iterates

% see VarAdaptAll for example of how to use

% timeType = 'eti' or 'ton'. If 'eti', estimates vs cycle time are calculated as well
% dt -- time step for calculation
% Ddt defines the width of the Gaussian prior
% logdat = 0 or 1 (shouldn't make a difference)
% output -- see lines 256-271. Set as [] to get everything outputted
% rescale_type: 'input' or 'output'. only guaranteed to work with 'input'
% multiType = 'add' or 'mult' sets the two-stimulus rate function to the mulitplication or linear combination


niter = 5; % number of times the alpha(t) estimate + fitting procedure is iterated
tolerance = 20; % if Delta(log-likelihood) drops below this number, the iteration procedure stops

% dalpha and alphamax set the alpha axis along which to find P(model|data)
dalpha = .5;
alphamax = 5;

if (isempty(var))
    btdstruct = BehaviorTriggeredData.prepVarianceSwitchingAnalysis_Gepner(btdstruct, opstruct);
    var = btdstruct.var;
    if(isfield(btdstruct, 'var_uv'))
        var_uv = btdstruct.var_uv;
    end
end

if(isempty(pd))
    pd_eti = BehaviorTriggeredData.createProblemStructForRateFunFitting (var, 'eti', opstruct.pdegree, opstruct.timeRange, [], []);
    pd_ton = BehaviorTriggeredData.createProblemStructForRateFunFitting (var, 'ton', opstruct.pdegree, opstruct.timeRange, [], []);
    pd_toff = BehaviorTriggeredData.createProblemStructForRateFunFitting (var, 'toff', opstruct.pdegree, opstruct.timeRange, [], []);
    
    if(strcmpi(timeType, 'eti'))
        pd = pd_eti;
    elseif(strcmpi(timeType, 'ton'))
        pd = pd_ton;
    elseif(strcmpi(timeType, 'toff'))
        pd = pd_toff;
    end

end

pd.maxreps = 2;
pd.pad = 0;


fn = fieldnames(pd);
for j = 1:length(fn)
    eval([fn{j} ' = pd.' fn{j} ';']);
end

if isempty(runEti) %#ok<*NODEF>
    runEti = runT;
end
if isempty(turnEti)
    turnEti = turnT;
end


if(strcmpi(timeType, 'ton'))
    separateExperiments = true;
    pd_ton = pd;
else
    separateExperiments = false;
end

if(strcmpi(timeType, 'eti') && (period+max(runT))<(max(runEti)))
    runTon = runT;
    turnTon = turnT;
    runT = runEti;
    turnT = turnEti;
end

tx = tx(:);

alpha = alpha_0;
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


if(exist('pd_ton'))
    runTon = pd_ton.runT;
    turnTon = pd_ton.turnT;
    turnTon = turnTon(tinds);
    runTon = runTon(rinds);
    [runTon, I] = sort(runTon, 'ascend');
    [turnTon, I] = sort(turnTon, 'ascend');
end
if(exist('runTon'))
    data.turnTon = turnTon;
    data.runTon = runTon;
end

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

data.tx = taxis;
data.txe = taxis;
data.deltaT = median(diff(var(1).fullensemble.eti));

    
funs.ratefun = ratefun;
funs.temporalratemod = [];
funs.params = params_0;
funs.tparams = tparams_0;
funs.tparams = [];
funs.gradlogratefun = gradlogratefun; %NxD gradient of log of rate
funs.hesslogratefun = hesslogratefun; %DxDxN hessian of log of rate;


if(size(alpha, 2)>2)
    alpha = alpha';
end
theta_0 = alpha(1,:);


if (separateExperiments)
    erange = [min(data.turnExpnum) max(data.turnExpnum)];
else
    erange = 1;
end

% set the rate function, and alpha_0 = 1
if(strcmpi(rescale_type, 'input'))
    rfun = @(alpha, rateP, dat) exp(polyval(rateP, alpha.*repmat(dat, size(alpha))));
    alphaT0 = ones(size(data.turnX));
    alphaR0 = ones(size(data.runX));
elseif(strcmpi(rescale_type, 'output'))
%     rfun = @(a, dat) a.*repmat( exp(polyval(funs.params_maxOutputRescaling, dat)), size(a));
    rfun = @(alpha, rateP, dat) alpha.*repmat( exp(polyval(rateP, dat)) - exp(rateP(end)), size(alpha)) + exp(rateP(end));
    alphaT0 = alphamax*ones(size(data.turnX));
    alphaR0 = alphamax*ones(size(data.runX));
end

% fit the rate function with alpha=1
r0 = log(length(pd.turnT)/(length(pd.runT)*data.deltaT));
[params{1}, nll(1)] = fitScaledRate (funs, data.turnX, alphaT0, data.runX, alphaR0, data.deltaT, r0, rescale_type, opstruct.pdegree, multiType);
    

% calculate alpha(t) and fit rate function with rescaled values, and iterate a few times to converge to the 'right' rate function
for i=2:(niter+1)
    scaled = 1;
    tic
    ratefun = @(alpha, dat) rfun(alpha, params{i-1}, dat);
    funs.params = params{i-1};
    
    
    if(strcmpi(timeType, 'ton'))
        tx = []; alpha = []; valpha = [];
        for j = min(erange):max(erange)
            alphaTon = BayesEstimateScalingFactor(ratefun, funs, data, theta_0, logdat, dalpha, alphamax, Ddt, j, timeType, opstruct.pdegree, multiType);
            alphaTon = normAlpha(alphaTon, data, []);
            tx = [tx alphaTon.tx];
            alpha = [alpha alphaTon.alpha];
            valpha = [valpha alphaTon.valpha];
        end
        [alphastruct.alpha_ton, alphastruct.valpha_ton, alphastruct.tx_ton] = meanEstinCycle(alpha, valpha, alphaTon.tx, tshift, period, opstruct.timeRange, timeType);
        alphastruct.alpha = alphastruct.alpha_ton;
        alphastruct.valpha = alphastruct.valpha_ton;
        alphastruct.tx = alphastruct.tx_ton;
        alphastruct = normAlpha(alphastruct, data, timeType);
        
        data = addStretchedValues(data, alphastruct, timeType);
        % fit rate function with new values of alpha(t)
        [params{i}, nll(i)] = fitStaticRate (funs, data, scaled, opstruct.pdegree, multiType);
        
    else
        
        tic
        alphastruct = BayesEstimateScalingFactor(ratefun, funs, data, theta_0, logdat, dalpha, alphamax, Ddt, [], timeType, opstruct.pdegree, multiType);
        toc
        alphastruct = normAlpha(alphastruct, data, timeType);
        
        data = addStretchedValues(data, alphastruct, timeType);
        % fit rate function with new values of alpha(t)
        [params{i}, nll(i)] = fitStaticRate (funs, data, scaled, opstruct.pdegree, multiType);
        
        [alphastruct.alpha_ton, alphastruct.valpha_ton, alphastruct.tx_ton] = meanEstinCycle(alphastruct.alpha, alphastruct.valpha, alphastruct.tx, tshift, period, opstruct.timeRange, timeType);
        alphastruct = normAlpha(alphastruct, data, 'ton');
    end

    alphastruct.params = params{i};
    alphastruct.nll = nll(i);
    alphastruct.Ddt = Ddt;
 
    fitstruct{i-1} = alphastruct;
    toc
    if( abs(nll(i)-nll(i-1)) < tolerance ) % stop iterating when log-likelihood converges
        break
    end
end

fitstruct{end}.scaledValues = ScaledData(fitstruct{end}, data, dt, timeType);
fitstruct{end}.nll_iterated = nll;
fitstruct{end}.params_iterated = params;


if(strcmpi(timeType, 'eti'))
    if(strcmpi(output, 'tx'))
        out = fitstruct{end}.tx;
    elseif(strcmpi(output, 'alpha'))
        out = fitstruct{end}.alpha;
    else
        out = fitstruct{end};
    end
else
    if(strcmpi(output, 'tx'))
        out = fitstruct{end}.tx_ton;
    elseif(strcmpi(output, 'alpha'))
        out = fitstruct{end}.alpha_ton;
    else
        out = fitstruct{end};
    end
end

end

function fitstruct = BayesEstimateScalingFactor(rfun, funs, data, alpha_0, logdat, dalpha, alphamax, Ddt, expnum, timeType, pdegree, multiType) % (ratefun, dlogratefun,ddlogratefun, turnT, turnX, runT, runX, deltaT, txe, Q_alpha, params, alpha_0, w_0)

existsAndDefault('alpha_0', 1);

rP = funs.params;

txe = data.txe;
theta = repmat(alpha_0, [length(txe), 1]);
dtheta = theta;
%theta = ones(size(txe));
innov = zeros(length(txe), 1);


if (existsAndDefault('expnum', []))
    tinds = data.turnExpnume == expnum;
    rinds = data.runExpnume == expnum;
else
    tinds = true(size(data.turnTe));
    rinds = true(size(data.runTe));
end

if(strcmpi(timeType, 'ton'))
    turnTe = data.turnTon(tinds);
    runTe = data.runTon(rinds);
else
    turnTe = data.turnTe(tinds);
    runTe = data.runTe(rinds);
end
turnEtie = data.turnTe(tinds);
turnXe = data.turnXe(tinds, :);
runEtie = data.runTe(rinds);
runXe = data.runXe(rinds, :);

nstim = size(runXe, 2);
% if(length(rP)>3 && nstim==2)
%     pdegree = 2;
% else
%     pdegree = 1;
% end

alphas = 0.05:dalpha:alphamax; %alpha axis for integration
[A1, A2] = meshgrid(alphas, alphas);

if(nstim==2)
    betamax = alphamax;
    beta_0 = alpha_0(2);
    alpha_0 = alpha_0(1);
    betas = 0.05:dalpha:betamax;
    [B1, B2] = meshgrid(betas, betas);
    [AB1, AB2] = meshgrid(alphas, betas);
    [BA1, BA2] = meshgrid(betas, alphas);
    if(length(Ddt)==3)
        sig2d = 2*[Ddt(1), Ddt(3); Ddt(3), Ddt(2)];
        if(det(sig2d)<=0)
            disp('sigma matrix must be positive definite');
            fitstruct.tx_ton = [];
            fitstruct.alpha_ton = [];
            fitstruct.isbad = 1;
            return
        end
    end
end


logPdat = @(a, tdat, rdat) sum( log(rfun(a,tdat)) )  - data.deltaT*sum(rfun(a,rdat)); %logP(data|model) = sum_turn log(rate) - dt * sum_noturn rate
Pdat = @(a, tdat, rdat) prod(data.deltaT*rfun(a,tdat), 1)  .* prod(exp(-data.deltaT*rfun(a,rdat)), 1); % (+ length(tdat)*log(dt)) P(data|model) 
Naa = @(a1,a2) normpdf(a1, a2, sqrt(2*Ddt)); % Diffusion prior on alpha P(alpha1 | alpha2), D=diffusion constant, dt = time step

if(nstim==2) % calculation for odor-light experiments
    
    scaledA = @(xd) alphas'.*repmat(xd(:,1)', length(alphas), 1);
    scaledB = @(xd) betas'.*repmat(xd(:,2)', length(betas), 1);
    scaledAmat = @(xd) repmat( reshape(scaledA(xd), length(alphas), 1, length(xd(:,1))), 1, length(betas), 1);
    scaledBmat = @(xd) repmat( reshape(scaledB(xd), [1, size(scaledB(xd))]), length(alphas), 1, 1);
    if(pdegree==1)
        Pdat2d = @(tdat, rdat) exp( sum( rP(1)*scaledAmat(tdat) + rP(2)*scaledBmat(tdat) + rP(3), 3) - data.deltaT*sum(exp( rP(1)*scaledAmat(rdat) + rP(2)*scaledBmat(rdat) + rP(3)), 3));
    elseif(pdegree==2)
        
            
        if(strncmpi(multiType, 'add', 3))
            Pdat2d = @(tdat, rdat) exp( sum( rP(2)* (cos(rP(1))*scaledAmat(tdat)+sin(rP(1))*scaledBmat(tdat)).^2 + rP(3)* (cos(rP(1))*scaledAmat(tdat)+sin(rP(1))*scaledBmat(tdat)) + rP(4), 3)...
                 - data.deltaT*sum(exp( rP(2)* (cos(rP(1))*scaledAmat(rdat)+sin(rP(1))*scaledBmat(rdat)).^2 + rP(3)* (cos(rP(1))*scaledAmat(rdat)+sin(rP(1))*scaledBmat(rdat)) + rP(4) ), 3) );   
        elseif(strncmpi(multiType, 'mult', 3))
            Pdat2d = @(tdat, rdat) exp( sum( rP(1)*scaledAmat(tdat).^2 + rP(2)*scaledAmat(tdat) + rP(3)*scaledBmat(tdat).^2 + rP(4)*scaledBmat(tdat) + rP(5), 3)...
             - data.deltaT*sum(exp( rP(1)*scaledAmat(rdat).^2 + rP(2)*scaledAmat(rdat) + rP(3)*scaledBmat(rdat).^2 + rP(4)*scaledBmat(rdat) + rP(5)), 3));

        end
        
    end
    Naa = @(a1, a2) normpdf(a1, a2, sqrt(2*Ddt(1))) ;
    Nbb = @(b1, b2) normpdf(b1, b2, sqrt(2*Ddt(2)));

    
end
    
%First time step

ti = findRangeInSortedData(turnTe, txe(1), txe(2));
ri = findRangeInSortedData(runTe, txe(1), txe(2));
td = turnXe(ti, :);
rd = runXe(ri, :);
xt = td.*repmat(theta(1, :), [size(td,1) 1]);
xr = rd.*repmat(theta(1, :), [size(rd,1) 1]);



if(nstim==1)
    if(logdat==1)
        Pax_0 = logPdat(alphas, xt, xr) + log(Naa(alphas, alpha_0));
        Pax_0 = Pax_0 - max(Pax_0); % extra term added to keep P from being too small
        omega = dalpha*trapz(exp(Pax_0));
        
        if(omega==0)
            Pax = exp(Pax_0);
        else
            Pax = exp(Pax_0)/omega;
        end
        Pa(1,:) = Pax;
    else
        Pax_0 = Pdat(alphas, xt, xr).* Naa(alphas, alpha_0);
        omega = dalpha*trapz(Pax_0);
        
        if(omega==0)
            Pax = Pax_0;
        else
            Pax = Pax_0/omega;
        end
        Pa(1,:) = Pax;
    end
elseif(nstim==2)
    
    if(length(Ddt)==2) % independent priors for a and b
        Pax_0 = Pdat2d(xt, xr).*Naa(alphas, alpha_0)'.*Nbb(betas, beta_0);
        
    elseif(length(Ddt)==3) %2D prior for (a,b)
        abmat = [AB1(:) AB2(:)];
        bamat = [BA1(:) BA2(:)];
        muab = [repmat(alpha_0, length(alphas)*length(betas), 1), repmat(beta_0, length(alphas)*length(betas), 1)];
        Pab_0 = mvnpdf(bamat, muab, sig2d);
        Pax_0 = Pdat2d(xt, xr).*reshape(Pab_0, size(BA1));
    end
    
    fitstruct.alphas = alphas;
    fitstruct.betas = betas;
    fitstruct.logPdat{1} = Pdat2d(xt, xr);
    omega = dalpha*dalpha*trapz(trapz(Pax_0));
    
    if(omega==0)
        Pax = Pax_0;
    else
        Pax = Pax_0/omega;
    end
    Pa{1} = Pax;
end


% fitstruct.Pdata(1,:) = Pdat(alphas, xt, xr)./omega;
% fitstruct.logPdata(1,:) = logPdat(alphas, xt, xr);
fitstruct.Pa(1,:) = Pa(1,:);

if(isempty(ri))
    inds = 1;
else
    inds = [];
end

for n = 2:(length(txe) - 1)
%    tEst = theta(j);
    
    ti = findRangeInSortedData(turnTe, txe(n), txe(n+1));
    ri = findRangeInSortedData(runTe, txe(n), txe(n+1));
    if (isempty(ri))
        inds = [inds n];
        theta(n+1,:) = theta(n,:);
        continue;
    end
    
    td = turnXe(ti, :);
    rd = runXe(ri, :);
    xt = td.*repmat(theta(n, :), [size(td,1) 1]);
    xr = rd.*repmat(theta(n, :), [size(rd,1) 1]);

    teti = turnEtie(ti);
    reti = runEtie(ri);
    
    if(nstim==1)
        
    Paa_n = dalpha*trapz( Naa(A1, A2).*repmat(Pax, length(alphas), 1), 2);
    
    if(logdat == 1)
       
        Pax_n = logPdat(alphas, xt, xr) + log(Paa_n');
        Pax_n = Pax_n - max(Pax_n);
        omega = dalpha*trapz( exp(Pax_n) );
        if(omega==0)
            Pax = exp(Pax_n);
        else
            Pax = exp(Pax_n)/omega;
        end
        Pa(n,:) = Pax;
    else
        Pax_n = Pdat(alphas, xt, xr).*Paa_n';
        omega = dalpha*trapz(Pax_n);
        if(omega==0)
            Pax = Pax_n;
        else
            Pax = Pax_n/omega;
        end
        Pa(n,:) = Pax;
    end
    fitstruct.Pa(n,:) = Pa(n,:);
    
    elseif(nstim==2)
        
        if(length(Ddt)==3) % 2-D prior
            
            rep = @(X) repelem(X(:), length(X(:)));
            rep2 = @(X) repmat(X(:), length(X(:)), 1);
            abmat = [rep2(AB1) rep2(AB2)];
            abmat2 = [rep(AB1) rep(AB2)];
            bamat = [rep2(BA1) rep2(BA2)];
            bamat2 = [rep(BA1) rep(BA2)];
            ab = mvnpdf(abmat, abmat2, sig2d);
            ba = mvnpdf(bamat, bamat2, sig2d);
            ab4 = reshape(ab, [size(AB1) size(AB2)]);
            ba4 = reshape(ba, [size(BA1) size(BA2)]);
            dat4 = repmat(Pax, [1, 1, size(Pax)]);
            Pab_n = dalpha*dalpha*trapz(trapz(ba4.*dat4, 2), 1);
            Pax_n = Pdat2d(xt, xr).*squeeze(Pab_n);
            omega = dalpha*dalpha*trapz(trapz(Pax_n));
            
        elseif(length(Ddt)==2) % independent 1-D priors
            
            % integrate alpha & beta together
            %         aa = Naa(A1, A2);
            %         bb = Nbb(B1, B2);
            %         aa2 = repmat(reshape(aa, length(alphas), 1, length(alphas), 1), 1, length(betas), 1, length(betas));
            %         bb2 = repmat(reshape(bb, 1, length(betas), 1, length(betas)), length(alphas), 1, length(alphas), 1);
            %         ab4 = aa2.*bb2;
            %         dat4 = repmat(Pax, [1, 1, size(Pax)]);
            %         Pab_n = dalpha*dalpha*trapz(trapz(ab4.*dat4, 2), 1);
            %         Pax_n = Pdat2d(xt, xr).*squeeze(Pab_n);
            
            % integrate alpha then beta
            Paa_n = squeeze(dalpha*trapz( repmat(Naa(A1,A2), 1, 1, length(betas)).*repmat( reshape(Pax,[1, size(Pax)]), length(alphas), 1, 1), 2));
            Pab_n = dalpha*trapz( repmat(reshape(Nbb(B1,B2), [1, size(Nbb(B1,B2))]), length(alphas),1, 1).*repmat( reshape(Paa_n,[size(Paa_n,1), 1, size(Paa_n,2)]), 1, length(betas), 1), 3);
            Pax_n = Pdat2d(xt, xr).*Pab_n;
            omega = dalpha*dalpha*trapz(trapz(Pax_n));
            
        end
        
        if(omega==0)
            Pax = Pax_n;
        else
            Pax = Pax_n/omega;
        end
        Pa{n} = Pax;
        fitstruct.Pa{n} = Pa{n};
        fitstruct.logPdat{n} = Pdat2d(xt, xr);
    end
    
% fitstruct.Pdata(n,:) = Pdat(alphas, xt, xr)./omega;
% fitstruct.logPdata(n,:) = logPdat(alphas, xt, xr);



end

m=1;
for n=1:(length(txe)-1)
    if(~ismember(n, inds))
        tEst(m) = txe(n);
        if(nstim==1)
%             if(logdat)
%                 Pan = exp(Pa(n,:))/(dalpha*trapz(exp(Pa(n,:))));
%                 meanEst(m) = trapz( alphas.* Pan ) / trapz(Pan);
%                 vEst(m) = trapz( alphas.^2.* Pa(n,:) ) / trapz(Pa(n,:)) - meanEst(m)^2;
%             else
                meanEst(m) = trapz( alphas.* Pa(m,:) ) / trapz(Pa(m,:));
                vEst(m) = trapz( (alphas-meanEst(m)).^2.*Pa(m,:) ) / trapz(Pa(m,:));
                vEst2(m) = trapz( alphas.^2.* Pa(n,:) ) / trapz(Pa(n,:)) - meanEst(m)^2;
%             end
        elseif(nstim==2)
            omega = trapz(trapz(Pa{n}));
            if(omega==0)
                meanEst(1, m) = dalpha*trapz( alphas'.*trapz(Pa{m}, 2) );
                meanEst(2, m) = dalpha*trapz( betas.*trapz(Pa{m}, 1) );
                vEst(1, 1, m) = dalpha*trapz( (alphas-meanEst(1,m)).^2.* Pa{m} );
                vEst(2, 2, m) = dalpha*trapz( (betas-meanEst(2,m))'.^2.* Pa{m} ) ;
                vEst(1,2,m) = trapz( (alphas-meanEst(1,m))'.*trapz((betas-meanEst(2,m)).*Pa{m}, 2) );
                vEst(2,1,m) = vEst(1,2,m);
            else
                meanEst(1,m) = trapz( alphas'.*trapz(Pa{m}, 2) )./omega;
                meanEst(2,m) = trapz( betas.*trapz(Pa{m}, 1) )./omega;
                vEst(1,1,m) = trapz( (alphas-meanEst(1,m))'.^2.* trapz(Pa{m},2) )./omega;
                vEst(2,2,m) = trapz( (betas-meanEst(2,m)).^2.* trapz(Pa{m},1) )./omega;
                vEst(1,2,m) = trapz( (alphas-meanEst(1,m))'.*trapz((betas-meanEst(2,m)).*Pa{m}, 2) )./omega;
                vEst(2,1,m) = vEst(1,2,m);
            end
        end
            
        m = m+1;
    end
end

fitstruct.Ddt = Ddt;
fitstruct.tx = tEst;
fitstruct.alpha = meanEst;
fitstruct.valpha = vEst;


end


function scaledValues = ScaledData(fitstruct, data, dt, timeType)

if(strcmpi(timeType, 'ton'))
    tx = fitstruct.tx_ton;
    alpha = fitstruct.alpha_ton;
else
    tx = fitstruct.tx;
    alpha = fitstruct.alpha;
end

for i=1:size(alpha,1)
    talpha(:,i) = interp1(tx, alpha(i,:), data.turnTe, 'spline');
    ralpha(:,i) = interp1(tx, alpha(i,:), data.runTe, 'spline');
end

scaledValues.deltaT = dt;
scaledValues.turnX = data.turnXe;
scaledValues.turnT = data.turnTe;
scaledValues.turnEti = data.turnEtie;
scaledValues.runX = data.runXe;
scaledValues.runT = data.runTe;
scaledValues.runEti = data.runEtie;

if(isfield(data, 'turnTon'))
scaledValues.turnTon = data.turnTon;
scaledValues.runTon = data.runTon;
end

scaledValues.alphaT = talpha;
scaledValues.alphaR = ralpha;
scaledValues.tval = data.turnXe .* talpha;
scaledValues.rval = data.runXe .* ralpha;

end

function [params, nLL, LL] = fitScaledRate(funs, tdata, talpha, rdata, ralpha, deltaT, r0, rescale_type, pdegree, multiType)

op = optimoptions('fminunc');
op.Display = 'off';
op.Algorithm = 'quasi-newton';

nstim = size(tdata, 2);


if(nstim==1)
    rate = @(p,xd) funs.ratefun(p, xd);
    if(pdegree==1)
        p0 = [0 r0];
    elseif(pdegree==2)
        p0 = [0 0 r0];
    end
elseif(nstim==2 && pdegree==2)
    if(strncmpi(multiType, 'add', 3))
        rate = @(p,xd) exp( p(2)* (cos(p(1))*xd(:,1) + sin(p(1))*xd(:,2)).^2 + p(3)* (cos(p(1))*xd(:,1) + sin(p(1))*xd(:,2)) + p(4) ) ;
        p0 = [pi/4 0 0 r0];
    elseif(strncmpi(multiType, 'mult', 3))
        rate = @(p,xd) exp( p(1)*xd(:,1).^2 + p(2)*xd(:,1) + p(3)*xd(:,2).^2 + p(4)*xd(:,2) + p(5) ) ;
        p0 = [0 0 0 0 r0];
    end
elseif(nstim==2 && pdegree==1)
    rate = @(p,xd) exp( p(1)*xd(:,1) + p(2)*xd(:,2) + p(3) );
    p0 = [0 0 r0];
end


if(strcmpi(rescale_type, 'input'))
    ratefun = @(alpha, rateP, dat) rate(rateP, alpha.*dat);
elseif(strcmpi(rescale_type, 'output'))
%     rfun = @(a, dat) a.*repmat( exp(polyval(funs.params_maxOutputRescaling, dat)), size(a));
    ratefun = @(alpha, rateP, dat) alpha.*( exp(polyval(rateP, dat)) - exp(rateP(end))) + exp(rateP(end));
end

nlogP = @(p) -sum(log(ratefun(talpha, p, tdata) ) ) + sum(ratefun(ralpha, p,rdata)*deltaT);

[params, nLL] = fminunc(nlogP, p0, op);
LL = -nLL + length(tdata)*log(deltaT);
end


function [params, ll] = fitStaticRate (funs, data, scaled, pdegree, multiType)
%function params = fitStaticRate (logratefun, params_0, turnT, turnX, runT, runX, tx, alpha)
%maximizes log(P(data|params)) = sum_turn logratefun(params, turnX*alpha) - sum_run exp(logratefun(params,runX*alpha)))*deltaT 

nstim = size(data.turnX, 2);

if(scaled)
    tval = data.tval;
    rval = data.rval;
else
    tval = data.turnXe;
    rval = data.runXe;
end
teti = data.turnEti;
reti = data.runEti;

op = optimoptions('fminunc');
op.Display = 'off';
op.Algorithm = 'quasi-newton';


if (isempty(funs.temporalratemod) || isempty(funs.tparams))
    
    funs.tparams = [];
    funs.temporalratemod = @(p,td) ones(size(td));
    
    if(nstim==1)
        ratefun = @(p,xd) funs.ratefun(p, xd);
    elseif(nstim==2 && pdegree==2)
        
        if(strncmpi(multiType, 'add', 3))
            ratefun = @(p,xd) exp( p(2)* (cos(p(1))*xd(:,1) + sin(p(1))*xd(:,2)).^2 + p(3)* (cos(p(1))*xd(:,1) + sin(p(1))*xd(:,2)) + p(4) ) ;
        elseif(strncmpi(multiType, 'mult', 3))
            ratefun = @(p,xd) exp( p(1)*xd(:,1).^2 + p(2)*xd(:,1) + p(3)*xd(:,2).^2 + p(4)*xd(:,2) + p(5) ) ;
        end
    elseif(nstim==2 && pdegree==1)
        ratefun = @(p,xd) exp( p(1)*xd(:,1) + p(2)*xd(:,2) + p(3) );
    end
    
%     if(length(funs.params)>5)
%         p0 = [funs.params(1:4) funs.params(end)];
%     end
        p0 = funs.params;
    
    nlogP = @(p) -sum(log(ratefun(p, tval)) ) + sum(ratefun(p,rval)*data.deltaT);
    [params, ll, exitflag] = fminunc(nlogP, p0, op);

    
else
    tparams = funs.tparams;
    params = funs.params;
    ratefun = @(p,tp,xd, td) funs.ratefun(p, xd) .* funs.temporalratemod(tp, td);
    nlogP = @(p)  -sum(log(ratefun(p, tparams, tval, teti)) ) + sum(ratefun(p,tparams, rval, reti)*data.deltaT);
    [params, ~, exitflag] = fminunc(nlogP, params, op);
    if (exitflag <= 0)
        warning ('static rate fit may not have converged!');
    end
    nlogP = @(tp)  -sum(log(ratefun(funs.params, tp, tval, teti)) ) + sum(ratefun(funs.params,tp, rval, reti)*data.deltaT);
    [tparams, ~, exitflag] = fminunc(nlogP, tparams, op);
    if (exitflag <= 0)
        warning ('static rate fit may not have converged!');
    end
    nlogP = @(p)  -sum(log(ratefun(p(1:length(params)), p((length(params)+1):end), tval, teti)) ) + sum(ratefun(p(1:length(params)), p((length(params)+1):end), rval, reti)*data.deltaT);
    p0 = [params tparams];
    [p, ll, exitflag] = fminunc(nlogP, p0, op);    
    params = p(1:length(funs.params));
    tparams = p((length(funs.params) + 1):end);
    
end



% ll = -ll + length(tval)*log(data.deltaT);
if (exitflag <= 0)
    warning ('static rate fit may not have converged!');
end
end


function data = addStretchedValues (data, fitstruct, timeType)

if(strcmpi(timeType, 'ton'))
    tx = fitstruct.tx_ton;
    alpha = fitstruct.alpha_ton;
else
    tx = fitstruct.tx;
    alpha = fitstruct.alpha;
end
    
if (length(fitstruct) == 1)
    [data.tval, data.rval] = stretchedValues (tx, alpha, data.turnT, data.turnX, data.runT, data.runX);
else
    data.tval = data.turnX;
    data.rval = data.runX;
    for j = 1:length(fitstruct)
        ti = data.turnExpnum == j;
        ri = data.runExpnum == j;
        [tval, rval] = stretchedValues(tx, alpha, data.turnT(ti), data.turnX(ti,:), data.runT(ri), data.runX(ri,:));
        data.tval(ti,:) = tval;
        data.rval(ri,:) = rval;
    end
end
end

function [tval, rval] = stretchedValues (tx, alpha, turnT, turnX, runT, runX)
for i=1:size(alpha,1)
    tval(:,i) = turnX(:,i) .* interp1(tx, alpha(i,:), turnT, 'spline');
    rval(:,i) = runX(:,i) .* interp1(tx, alpha(i,:), runT, 'spline');
end
end



function fitstruct = normAlpha (fitstruct, data, timeType, expnum)
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
    alpha(isnan(alpha)) = 0;           %jwolk
    norm_factor = sum(alpha.*repmat(wt, [1 size(alpha,2)]),1);
    alpha = alpha ./ repmat(norm_factor, [size(alpha,1) 1]);
    fitstruct.alpha = alpha';
end

end
