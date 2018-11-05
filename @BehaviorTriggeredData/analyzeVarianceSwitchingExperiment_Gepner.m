function btdstruct = analyzeVarianceSwitchingExperiment_Gepner(btdstruct, opstruct, qfactor)
%btdstruct = analyzeVarianceSwitchingExperiment(btdstruct, opstruct)
%

if (nargin < 1)
    opstruct.stim(1).gqname = 'led1ValDiff';
    opstruct.stim(1).iswn = true;
    opstruct.stim(1).isstep = false;
    opstruct.stim(1).ramptype = 'square';
    opstruct.stim(1).switch = true;
    opstruct.stim(1).period = [];
    opstruct.stim(2).gqname = 'led2ValDiff';
    opstruct.stim(2).iswn = true;
    opstruct.stim(2).isstep = false;
    opstruct.stim(2).ramptype = 'square';
    opstruct.stim(2).switch = false;
    opstruct.stim(2).period = [];
    opstruct.timeRange = [60 Inf];
    opstruct.numTimeBins = 120; %time bin size = period/numTimeBins
    opstruct.redogq = false;
    opstruct.numLxBins = 10;
    opstruct.kernelTime = 7;
    opstruct.kernelDt = 0.1;
    opstruct.adaptationTime = 10; %first T seconds to exclude from high/low cycles in analysis
    opstruct.Q_alpha = 0.2;
    opstruct.pdegree = 1;
    btdstruct = opstruct;
    disp ('ramp types are square, triangle, sine, and constant');
    return;
end

if (nargin < 2)
    error ('first argument must be a btdstruct and second argument contains options');
end
% 
% opstruct.kalmanNoiseDiffusion1D = qfactor*eye(3);
% opstruct.kalmanNoiseDiffusion2D = qfactor*eye(6);

if(~isfield(btdstruct,'varops'))
    btdstruct.varops = opstruct;
else
    btdstruct.varops.kalmanNoiseDiffusion1D = qfactor*eye(3);
    btdstruct.varops.kalmanNoiseDiffusion2D = qfactor*eye(6);
    opstruct = btdstruct.varops;
end


n = find([opstruct.stim.varswitch]);
gq_switch = opstruct.stim(n).gqname;

% if(~isfield(btdstruct, 'btd'))
%     
%     period = btdstruct.varops.stim.period;
%     
%     turnToff = btdstruct.var.turn.toff;
%     turnTon = btdstruct.var.turn.ton;
%     turnX = btdstruct.var.turn.x_conv;
%     runToff = btdstruct.var.noturn.toff;
%     runTon = btdstruct.var.noturn.ton;
%     runX = btdstruct.var.noturn.x_conv;
%     samplingInterval = median(diff(btdstruct.var.fullensemble.eti));
%     qrate = btdstruct.varops.kalmanNoiseDiffusion;
%     
%     btdstruct.var.rateEst_on = estimateRateVsTime(turnTon, turnX, runTon, runX, period, samplingInterval, qrate);
%     btdstruct.var.rateEst_off = estimateRateVsTime(turnToff, turnX, runToff, runX, period, samplingInterval, qrate);
%     return;
% end


if (~opstruct.stim(1).iswn || strcmpi(opstruct.stim(1).ramptype, 'constant'))
    error ('reorder options so that first argument is a white noise signal with time varying variance');
end


for j = find([opstruct.stim.iswn])
    gqname = opstruct.stim(j).gqname;
    
    if(opstruct.stim(j).isstep) %take out the step from the stimulus
        gq_nostep = strcat(gqname, '_nostep');
        for i=1:length(btd)
            ind = btd(i).findField(gqname);
            old = btd(i).glt(ind);
            new = old;
            new.fieldname = gq_nostep;
            new.yData(abs(new.yData)>7*std(new.yData))=0;
            btdstruct.btd(i).glt(end+1) = new;
        end
        gqname = gq_nostep;
    end
    
    if (opstruct.stim(j).varswitch == 1)
        if (any(btdstruct.btd.findField([gqname '_var_high']) <= 0))
            btdstruct.btd = btdstruct.btd.addVarianceGQs(gqname, opstruct.stim(j).ramptype);
        end
    end
    if ( j<3 && (opstruct.redogq || any(btdstruct.btd.findField([gqname 'TurnLin']) <= 0)))
        btdstruct = addScaledVarianceTurnLin (btdstruct, opstruct, j);
    end
    
end


% for j = find(~[opstruct.stim.iswn])
%     if (opstruct.redogq || any(btdstruct.btd.findField([opstruct.stim(j).gqname '_ton']) <= 0))
%         btdstruct.btd = btdstruct.btd.addTonToffGQs(opstruct.stim(j).gqname, opstruct.stim(j).ramptype, 'fixedPeriod', opstruct.stim(j).period);
%     end
% end
% 
% if (any(btdstruct.btd.findField([gq_switch '_ratePredROG_var_low']) <= 0) || any(btdstruct.btd.findField([gq_switch '_ratePredROG_var_low_last'])<=0) )
%     btdstruct = ratefunGQVariance(btdstruct, opstruct);
% end


if (length(opstruct.stim) == 1)
    btdstruct = analyzeVariance1Source (btdstruct);
%     btdstruct = VarInnovations(btdstruct);
    return
end


if (nnz([opstruct.stim.iswn]) > 1)
    btdstruct = analyzeVariancePlusNoise (btdstruct, n);
%     btdstruct = VarInnovations(btdstruct);
    return;
end
end

function btdstruct = addScaledVarianceTurnLin (btdstruct,opstruct,j)

btd = btdstruct.btd;
gqname = opstruct.stim(j).gqname;
if(opstruct.stim(j).isstep)
    gqname = [gqname '_nostep'];
end
kname = [gqname 'TurnLin'];
kname_low = [gqname 'TurnLin_var_low'];
kname_high = [gqname 'TurnLin_var_high'];

kernelTime = opstruct.kernelTime;
kernelDt = opstruct.kernelDt;
btdstruct.kernelTime = kernelTime;
btdstruct.kernelDt = kernelDt;

taxis = -(0:kernelDt:kernelTime);
    
if( nnz([opstruct.stim.iswn]) == 1 && opstruct.stim(j).varswitch == 1 && strcmpi(opstruct.stim(j).ramptype, 'square'))
    
    
    stimTon = btd.behaviorTriggeredDataMatrix ('all', '', [gqname '_var_ton'],0);
    stimToff = btd.behaviorTriggeredDataMatrix ('all', '', [gqname '_var_toff'],0);
    turnTon = btd.behaviorTriggeredDataMatrix ('turn', 'start', [gqname '_var_ton'],0);
    turnToff = btd.behaviorTriggeredDataMatrix ('turn', 'start', [gqname '_var_toff'],0);
    
    islow = turnToff < turnTon & turnToff > 0;
    ishigh = turnToff > turnTon & turnToff > 0;
    islow_all = stimToff < stimTon & stimToff > 0;
    ishigh_all = stimToff > stimTon & stimToff > 0;
    
    %using only one kernel

    btdstruct.linname{j} = kname;
    [convkernel, btd] = btdstruct.btd.createBTAKernel('turn', 'start', gqname, kernelTime, kernelDt, 'newFieldName',kname, 'abbott', true);
    stimInput = btd.behaviorTriggeredDataMatrix ('all', '', kname, 0);
    turnInput = btd.behaviorTriggeredDataMatrix('turn', 'start', gqname, taxis);
    
    inds_low = islow_all & stimToff > kernelTime; %guarantee only in low variance portion
    inds_high = ishigh_all & stimTon>kernelTime;
    sl = std(stimInput(inds_low));
    sh = std(stimInput(inds_high));
    s = sl;
    
    btdstruct.kernelScaling(j) = s;
    convkernel = convkernel / s;
    
    
    % using high and low kernels
    
%     aa = [btd.all];
%     alleti = [aa.eti];
%     tt = [btd.turn]';
%     turneti = [tt.start_eti]';
%     
%     mintime = min(alleti(stimTon > 0 & stimToff > 0));
%     turnvalid = turneti >= mintime;
%     allvalid = alleti >= mintime;
% 
%     ishigh = ishigh(turnvalid);
%     islow = islow(turnvalid);
%     
%     dm_all = turnInput(turnvalid, :);
%     
%     dm_high = dm_all(ishigh, :);
%     dm2 = dm_high; dm2(~isfinite(dm2)) = 0;
%     bta_high = sum(dm2)./sum(isfinite(dm_high));
%     
%     dm_low = dm_all(islow, :);
%     dm2 = dm_low; dm2(~isfinite(dm2)) = 0;
%     bta_low = sum(dm2)./sum(isfinite(dm_low));
%     bta_low = bta_low * max(abs(bta_high))/max(abs(bta_low));
%     
%     abbottfun = @(x,t) (x(3)*exp(-x(2)*t)/(x(1)-x(2)) - x(5)*exp(-x(4)*t)/(x(1)-x(4)) + (x(3)*(x(4)-x(1)) - x(5)*(x(2)-x(1)))*exp(-x(1)*t)/((x(1)-x(2))*(x(1)-x(4)))) .* (t <= 0);
%     
%     [~,I] = max(abs(bta_low));
%     mv = bta_low(I);
%     x0 = [-1.5 -1 mv -0.5 mv*1.5];
%     op = optimset('lsqcurvefit');
%     op.MaxFunEvals = 1E4;
%     op.MaxIter = 1e4;
%     x = lsqcurvefit(abbottfun, x0, taxis, bta_low,[],[],op);
%     convkernel_low = abbottfun(x,taxis);
%     
%     
%     [~,I] = max(abs(bta_high));
%     mv = bta_high(I);
%     x0 = [-1.5 -1 mv -0.5 mv*1.5];
%     op = optimset('lsqcurvefit');
%     op.MaxFunEvals = 1E5;
%     op.MaxIter = 1e4;
%     x = lsqcurvefit(abbottfun, x0, taxis, bta_high,[],[],op);
%     convkernel_high = abbottfun(x,taxis);
%     
%     convkernel_low = convkernel_low * max(abs(convkernel))/max(abs(convkernel_low));
%     convkernel_high = convkernel_high * max(abs(convkernel))/max(abs(convkernel_high));
%     
%     btdstruct.btd = btdstruct.btd.addConvolvedFields(gqname, kname_low, convkernel_low, kernelDt, 'scaleToSqr', true);
%     btdstruct.btd = btdstruct.btd.addConvolvedFields(gqname, kname_high, convkernel_high, kernelDt, 'scaleToSqr', true);
%     
%     btdstruct.convkernels_low{j} = convkernel_low;
%     btdstruct.convkernels_high{j} = convkernel_high;
    
else
    btd = btdstruct.btd;
    btdstruct.kernelScaling(j) = opstruct.kernelScaling(j);
    convkernel = opstruct.convkernels{j};
end

btdstruct.btd = btdstruct.btd.addConvolvedFields(gqname, kname, convkernel, kernelDt, 'scaleToSqr', true);

btdstruct.convkernels{j} = convkernel;


end


function rateEstimate = estimateRateVsTime (turnT, turnT_modT, turnX, runT, runT_modT, runX, allT, allT_modT, allX, period, samplingInterval, qrate, opstruct)

turnT = turnT';
turnX = turnX';
allT = allT';
runT = runT';
runX = runX';

allX = [turnX, runX];


ratefun_names = {'rateExp1', 'rateExp2', 'rateExp1_rescale', 'rateExp2_rescale'};

if(length(qrate)==3)
    nstim = 1;
    rateExp1 = @(x, xdata) exp( x(1)*xdata + x(2) );
    rateExp2 = @(x, xdata) exp( x(1)*xdata.^2 + x(2)*xdata + x(3) );
elseif(length(qrate)==6)
    nstim = 2;
    rateExp1 = @(x, xdata) exp(x(1) * xdata(1,:) + x(2) * xdata(2,:) + x(3));
    rateExp2 = @(x, xdata) exp( x(1)*xdata(1,:).^2 + x(2)*xdata(2,:).^2 + x(3)*xdata(1,:).*xdata(2,:) + x(4)*xdata(1,:) + x(5)*xdata(2,:) + x(6) );
end

ratefuns = {rateExp1, rateExp2};

% using exponential(k==1) and ROG (k==2) fits for r(x)
% for k=1:2
%
%
%     [~, ~, theta0,w0] = PointProcessFilter_Gepner(ratefun_names{k}, [], turnT, turnX, runT, runX, qrate, samplingInterval,(0:1:period));
%     [~, ~, theta0,w0,etx] = PointProcessFilter_Gepner(ratefun_names{k}, [], turnT, turnX, runT, runX, qrate, samplingInterval,(0:.5:period),theta0(:,end), w0(:,:,end));
%     ind0 = find(etx >= 0.75*period, 1, 'first');
%     [dkl, dkl_eb, theta, w, est_tx, diagnostics] = PointProcessFilter_Gepner(ratefun_names{k}, [], turnT, turnX, runT, runX, qrate, samplingInterval,0:dtfinal:period,theta0(:,ind0), w0(:,:,ind0));
%     [theta_s, w_s] = applySmoothingGaussianFilters(theta, w, qrate*dtfinal);
%
%
%
%     inds = est_tx >= 0 & est_tx <= period;
%     rateEstimate.rateExp{k}.tx = est_tx(inds);
%     rateEstimate.rateExp{k}.theta = theta(:,inds);
%     rateEstimate.rateExp{k}.theta_s = theta_s(:,inds);
%     rateEstimate.rateExp{k}.w = w(:,:,inds);
%     rateEstimate.rateExp{k}.w_s = w_s(:,:,inds);
%
%     if (size(theta,1)==3)
%         rateEstimate.rateExp{k}.dkl = dkl(inds);
%         rateEstimate.rateExp{k}.dkl_eb = dkl_eb(inds);
% %         rateEstimate.rateExp{k}.dkl_s = dkl_s(inds);
% %         rateEstimate.rateExp{k}.dkl_eb_s = dkl_eb_s(inds);
%
%     elseif (size(theta,1)==6)
%
%         rateEstimate.rateExp{k}.dkl = dkl{1}(inds); rateEstimate.dkl_s = dkl_s{1}(inds);
%         rateEstimate.rateExp{k}.dkl_l = dkl{2}(inds); rateEstimate.dkl_l_s = dkl_s{2}(inds);
%         rateEstimate.rateExp{k}.dkl_o = dkl{3}(inds); rateEstimate.dkl_o_s = dkl_s{3}(inds);
%
%         rateEstimate.rateExp{k}.dkl_eb = dkl_eb(inds);
%
%
%     end
%
%     rateEstimate.rateExp{k}.logP = diagnostics.logP(inds);
%     rateEstimate.rateExp{k}.innovation = diagnostics.innovation(inds);
%     rateEstimate.rateExp{k}.innovation_scaled = diagnostics.innovation_scaled(inds);
%     rateEstimate.rateExp{k}.innovation_scaled_alt = diagnostics.innovation_scaled_alt(inds);
% end


% Gain rescaling only: r(t) = r(a(t)*x) or r(t) = r(a(t)*x1, b(t)*x2)

tolerance = 1e-4;

dtfinal = 0.2;

r0 = log(size(turnX, 2))-log(samplingInterval*size(runX, 2));

tmin = min(allT);
tmax = max(allT);

tx1 = (tmin:1:tmax);
tx2 = (tmin:0.5:tmax);
txf = tmin:dtfinal:tmax;

taxis1 = 0:1:period;
taxis2 = 0:.5:period;
taxisf = 0:dtfinal:period;

histaxis = linspace(-5, 5, 20);

for k=1:1
    
    clear params fitParams dParams dt alpha beta theta0 theta;
    
    %initialize rate-function fit-parameters
    if (nstim==1 && k==1)
        p0 = [0 r0 0];
        fitParams = fitRate(ratefuns{k}, p0(1:2), turnX, runX, samplingInterval, opstruct.numLxBins );
        p1 = [fitParams 0];
    elseif(nstim==1 && k==2)
        p0 = [0 0 r0];
        fitParams = fitRate(ratefuns{k}, p0, turnX, runX, samplingInterval, opstruct.numLxBins );
        p1 = fitParams;
    elseif(nstim==2 && k==1)
        p0 = [0 0 r0 0 0 0];
        fitParams = fitRate(ratefuns{k}, p0(1:3), turnX, runX, samplingInterval, opstruct.numLxBins );
        p1 = [fitParams 0 0 0];
    elseif(nstim==2 && k==2)
        p0 = [0 0 0 0 0 r0];
        fitParams = fitRate(ratefuns{k}, p0, turnX, runX, samplingInterval, opstruct.numLxBins );
        p1 = fitParams;
    end
    
    if(max(allT) == max(allT_modT))
        
        turnX_rescaled = zeros(size(turnX)) ;
        runX_rescaled = zeros(size(runX)) ;
        allX_rescaled = zeros(size(allX)) ;

        params{1} = p1;
        dParams{1} = p1 - p0;
        
        %estimate time-evolution of rescaling parameters:
        %(theta(1) for 1-stim. expts., theta(1) and theta(2) for 2-stim. expts.
        % then fit rate-function params again using rescaled stimulus

        
        n=1;
        while( sum(dParams{n}.^2) > tolerance )
            clear dt alpha beta theta0 theta w0 w turnX_rescaled runX_rescaled allX_rescaled;
            
            tx = tx1;
            taxis = taxis1;
            [~, ~, theta] = PointProcessFilter_Gepner(ratefun_names{k+2}, p1', turnT, turnX, runT, runX, qrate, samplingInterval,tx, p0');
            
            
            alpha = theta(1,:)./mean(theta(1,:));
            %[~, malpha] = meanyvsx(mod(tx, period), alpha, binEdgesFromCenters(taxis));
            
            if(nstim==2)
                beta = theta(2,:)./mean(theta(2,:));
                %[~, mbeta] = meanyvsx(mod(tx, period), beta, binEdgesFromCenters(taxis));
            end
            
            [~,~,turnBin] = histcounts(turnT_modT,  binEdgesFromCenters(tx));
            [~,~,runBin] = histcounts(runT_modT,  binEdgesFromCenters(tx));
            [~,~,allBin] = histcounts(allT_modT,  binEdgesFromCenters(tx));
            
            
            dt = diff(binEdgesFromCenters(tx));
            for j = 1:length(dt)
                
                turninds = turnBin == j;
                runinds = runBin == j;
                allinds = allBin == j;
                turnX_rescaled(1,turninds) = alpha(j) * turnX(1,turninds) ;
                runX_rescaled(1,runinds) = alpha(j) * runX(1,runinds);
                if(nstim==2)
                    turnX_rescaled(2,turninds) = beta(j) * turnX(2,turninds);
                    runX_rescaled(2,runinds) = beta(j) * runX(2,runinds);
                end
            end
            
            if(k==1)
                fitParams = fitRate(ratefuns{k}, p1(1:nstim+1), turnX_rescaled, runX_rescaled, samplingInterval, opstruct.numLxBins);
                if(nstim==1), p2 = [fitParams 0];
                elseif(nstim==2), p2 = [fitParams 0 0 0];
                end
            elseif(k==2)
                fitParams = fitRate(ratefuns{k}, p1, turnX_rescaled, runX_rescaled, samplingInterval, opstruct.numLxBins);
                p2 = fitParams;
            end
            
            dParams{n+1} = p2 - p1;
            params{n+1} = p2;
            
            p0 = theta(:,end);
            p1 = p2;
            n = n+1;
            
        end
        rateEstimate.rateExp{k}.theta0 = theta;
        rateEstimate.rateExp{k}.rateParams = params;
        
    else
        
        
        rP = opstruct.rateExp{k}.rateParams;
        
        tx = tx2;
        [~, ~, theta, w, ~, diagnostics] = PointProcessFilter_Gepner(ratefun_names{k+2}, rP', turnT, turnX, runT, runX, qrate, samplingInterval,tx, p0');
        
%         t0 = theta(:,end);
%         
%         tx = txf;
%         [~, ~, theta, w, ~, diagnostics] = PointProcessFilter_Gepner(ratefun_names{k+2}, rP', turnT, turnX, runT, runX, qrate, samplingInterval,tx, t0');
%         
        [theta_s, w_s] = applySmoothingGaussianFilters(theta, w, qrate*median(diff(tx)));
        
        
        tshift = opstruct.timeRange(1) - period/2;
        if(tshift<0)
            tshift  = -tshift;
        end
        
        [thetamodT, thetamodT_eb, tmodT] = etiThetaToCycleTheta(theta', w, tx, -tshift, period);
        [thetamodT_s, thetamodT_eb_s, ~] = etiThetaToCycleTheta(theta_s', w_s, tx, -tshift, period);
        
        rateEstimate.rateExp{k}.tx = tx;
        rateEstimate.rateExp{k}.theta = theta;
        rateEstimate.rateExp{k}.w = w;
        rateEstimate.rateExp{k}.theta_s = theta_s;
        rateEstimate.rateExp{k}.w_s = w_s;
        
        rateEstimate.rateExp{k}.tmodT = tmodT;
        rateEstimate.rateExp{k}.alpha = thetamodT(:,1)'./mean(thetamodT(:,1));
        rateEstimate.rateExp{k}.alpha_eb = squeeze(thetamodT_eb(1,1,:))';
        rateEstimate.rateExp{k}.alpha_s = thetamodT_s(:,1)'./mean(thetamodT_s(:,1));
        rateEstimate.rateExp{k}.alpha_eb_s = squeeze(thetamodT_eb_s(1,1,:))';
        
        rateEstimate.rateExp{k}.logP_rescaling = diagnostics.logP;
        rateEstimate.rateExp{k}.innovations = diagnostics.innovation;
        rateEstimate.rateExp{k}.innovations_scaled = diagnostics.innovation_scaled;
        
        if(nstim==2)
            rateEstimate.rateExp{k}.beta = thetamodT(:,2)'./mean(thetamodT(:,2));
            rateEstimate.rateExp{k}.beta_eb = squeeze(thetamodT_eb(2,2,:))';
            rateEstimate.rateExp{k}.beta_s = thetamodT_s(:,2)'./mean(thetamodT_s(:,2));
            rateEstimate.rateExp{k}.beta_eb_s = squeeze(thetamodT_eb_s(2,2,:))';
        end
    end
end

end


function btdstruct = analyzeVariance1Source (btdstruct, whichStim)

existsAndDefault ('whichStim', 1);

btd = btdstruct.btd;
gq = btdstruct.varops.stim(whichStim).gqname;
gq_conv = strcat(gq, 'TurnLin');

if(strcmpi(btdstruct.varops.stim(whichStim).ramptype, 'square'))
    gq_varh = strcat(gq, '_var_high');
    gq_varl = strcat(gq, '_var_low');
else
    gq_varh = strcat(gq, '_var_rising');
    gq_varl = strcat(gq, '_var_falling');
end

gq_ton = [gq,'_var_ton'];
gq_toff = [gq,'_var_toff'];


if (isempty(btdstruct.varops.stim(whichStim).period))
    btdstruct.varops.stim(whichStim).period = max(btd.behaviorTriggeredDataMatrix('all',[], gq_ton, 0));
end
period = btdstruct.varops.stim(whichStim).period;

if (~isfield(btdstruct, 'var') || length(btdstruct.var) < whichStim)
    
    aa = [btd.all]';
    alleti = [aa.eti]';
    tt = [btd.turn]';
    turneti = [tt.start_eti]';
    isrun = [aa.isrun]';
    
    all_ton = btd.behaviorTriggeredDataMatrix('all', 'start', gq_ton, 0);
    all_toff = btd.behaviorTriggeredDataMatrix('all', 'start', gq_toff, 0);
    turn_ton = btd.behaviorTriggeredDataMatrix('turn', 'start', gq_ton, 0);
    turn_toff = btd.behaviorTriggeredDataMatrix('turn', 'start', gq_toff, 0);
    
    mintime = min(alleti(all_ton > 0 & all_toff > 0));
    btdstruct.varops.timeRange(1) = max(btdstruct.varops.timeRange(1), mintime);
    
    turnvalid = turneti >= min(btdstruct.varops.timeRange) & turneti < max(btdstruct.varops.timeRange);
    allvalid = alleti >= min(btdstruct.varops.timeRange) & alleti < max(btdstruct.varops.timeRange);
    
    turndata = btd.behaviorTriggeredDataMatrix('turn', 'start', gq_conv, 0);
    alldata =  btd.behaviorTriggeredDataMatrix('all', 'start', gq_conv, 0);

    btdstruct.var(whichStim).all.x_conv = alldata(allvalid);
    btdstruct.var(whichStim).all.eti = alleti (allvalid);
    btdstruct.var(whichStim).all.ton = all_ton (allvalid);
    btdstruct.var(whichStim).all.toff = all_toff (allvalid);
    
    
    btdstruct.var(whichStim).run.x_conv = alldata(isrun & allvalid);
    btdstruct.var(whichStim).run.eti = alleti (isrun & allvalid);
    btdstruct.var(whichStim).run.ton = all_ton (isrun & allvalid);
    btdstruct.var(whichStim).run.toff = all_toff (isrun & allvalid);
    
    btdstruct.var(whichStim).turn.x_conv = turndata(turnvalid);
    btdstruct.var(whichStim).turn.eti = turneti (turnvalid);
    btdstruct.var(whichStim).turn.ton = turn_ton (turnvalid);
    btdstruct.var(whichStim).turn.toff = turn_toff (turnvalid);
    
%     if(strcmpi(btdstruct.varops.stim(whichStim).ramptype, 'square'))
%         
%         turndata_low = btd.behaviorTriggeredDataMatrix('turn', 'start', [gq_conv '_var_low'], 0);
%         turndata_high = btd.behaviorTriggeredDataMatrix('turn', 'start', [gq_conv '_var_high'], 0);
%         alldata_low =  btd.behaviorTriggeredDataMatrix('all', 'start', [gq_conv '_var_low'], 0);
%         alldata_high =  btd.behaviorTriggeredDataMatrix('all', 'start', [gq_conv '_var_high'], 0);
%         btdstruct.var(whichStim).all.x_conv_low = alldata_low(allvalid);
%         btdstruct.var(whichStim).all.x_conv_high = alldata_high(allvalid);
%         btdstruct.var(whichStim).run.x_conv_low = alldata_low(isrun & allvalid);
%         btdstruct.var(whichStim).run.x_conv_high = alldata_high(isrun & allvalid);
%         btdstruct.var(whichStim).turn.x_conv_low = turndata_low(turnvalid);
%         btdstruct.var(whichStim).turn.x_conv_high = turndata_high(turnvalid);
%     end

end


turnToff = btdstruct.var(whichStim).turn.toff;
turnTon = btdstruct.var(whichStim).turn.ton;
turnEti = btdstruct.var(whichStim).turn.eti;
turnX = btdstruct.var(whichStim).turn.x_conv;
runToff = btdstruct.var(whichStim).run.toff;
runTon = btdstruct.var(whichStim).run.ton;
runEti = btdstruct.var(whichStim).run.eti;
runX = btdstruct.var(whichStim).run.x_conv;
allToff = btdstruct.var(whichStim).all.toff;
allTon = btdstruct.var(whichStim).all.ton;
allEti = btdstruct.var(whichStim).all.eti;
allX = btdstruct.var(whichStim).all.x_conv;
samplingInterval = median(diff(btdstruct.var(whichStim).all.eti));
qrate = btdstruct.varops.kalmanNoiseDiffusion1D;

if(~isfield(btdstruct.varops, 'rateExp'))
    rateEst_off = estimateRateVsTime(turnToff, turnToff, turnX, runToff, runToff, runX, allToff, allToff, allX, period, samplingInterval, qrate, btdstruct.varops);
    for k=1:length(rateEst_off.rateExp)
        btdstruct.varops.rateExp{k}.rateParams = rateEst_off.rateExp{k}.rateParams{end};
        btdstruct.varops.rateExp{k}.theta0 = rateEst_off.rateExp{k}.theta0;
    end
end

btdstruct.var(whichStim).rateEst_off = estimateRateVsTime(turnEti, turnToff, turnX, runEti, runToff, runX, allEti, allToff, allX, period, samplingInterval, qrate, btdstruct.varops);
for k=1:length(btdstruct.var(whichStim).rateEst_off.rateExp)
    btdstruct.var(whichStim).rateEst_off.rateExp{k}.rateParams = btdstruct.varops.rateExp{k}.rateParams;
    btdstruct.var(whichStim).rateEst_off.rateExp{k}.theta0 = btdstruct.varops.rateExp{k}.theta0;
end


%btdstruct.var(whichStim).rateEst_on = estimateRateVsTime(turnTon, turnX, runTon, runX, allTon, allX, period, samplingInterval, qrate, btdstruct.varops);

% if(strcmpi(btdstruct.varops.stim(whichStim).ramptype, 'square'))
%     runX_low = btdstruct.var(whichStim).noturn.x_conv_low;
%     runX_high = btdstruct.var(whichStim).noturn.x_conv_high;
%     turnX_low = btdstruct.var(whichStim).turn.x_conv_low;
%     turnX_high = btdstruct.var(whichStim).turn.x_conv_high;
%     btdstruct.var(whichStim).rateEst_on_lowK = estimateRateVsTime(turnTon, turnX_low, runTon, runX_low, period, samplingInterval, qrate, opstruct);
%     btdstruct.var(whichStim).rateEst_off_lowK = estimateRateVsTime(turnToff, turnX_low, runToff, runX_low, period, samplingInterval, qrate, opstruct);
%     
%     btdstruct.var(whichStim).rateEst_on_highK = estimateRateVsTime(turnTon, turnX_high, runTon, runX_high, period, samplingInterval, qrate, opstruct);
%     btdstruct.var(whichStim).rateEst_off_highK = estimateRateVsTime(turnToff, turnX_high, runToff, runX_high, period, samplingInterval, qrate, opstruct);
% end


end




function btdstruct = analyzeVariancePlusNoise (btdstruct, whichStim)

btd = btdstruct.btd;    

kernelTime = btdstruct.varops.kernelTime;

gq_switch = btdstruct.varops.stim(whichStim).gqname;

if(strcmpi(btdstruct.varops.stim(whichStim).ramptype, 'square'))
    gq_varh = strcat(gq_switch, '_var_high');
    gq_varl = strcat(gq_switch, '_var_low');
else
    gq_varh = strcat(gq_switch, '_var_rising');
    gq_varl = strcat(gq_switch, '_var_falling');
end

gq_ton = [gq_switch,'_var_ton'];
gq_toff = [gq_switch,'_var_toff'];

aa = [btd.all]';
alleti = [aa.eti]';
tt = [btd.turn]';
turneti = [tt.start_eti]';
isrun = [aa.isrun]';

all_ton = btd.behaviorTriggeredDataMatrix('all', 'start', gq_ton, 0);
all_toff = btd.behaviorTriggeredDataMatrix('all', 'start', gq_toff, 0);
turn_ton = btd.behaviorTriggeredDataMatrix('turn', 'start', gq_ton, 0);
turn_toff = btd.behaviorTriggeredDataMatrix('turn', 'start', gq_toff, 0);

mintime = min(alleti(all_ton > 0 & all_toff > 0));
btdstruct.varops.timeRange(1) = max(btdstruct.varops.timeRange(1), mintime);

turnvalid = turneti >= min(btdstruct.varops.timeRange) & turneti < max(btdstruct.varops.timeRange);
allvalid = alleti >= min(btdstruct.varops.timeRange) & alleti < max(btdstruct.varops.timeRange);

if (isempty(btdstruct.varops.stim(whichStim).period))
    btdstruct.varops.stim(whichStim).period = max(btd.behaviorTriggeredDataMatrix('all',[], gq_ton, 0));
end

period = btdstruct.varops.stim(whichStim).period;


stimnames = {btdstruct.varops.stim(1).gqname, btdstruct.varops.stim(2).gqname};
kname1 = [stimnames{1} 'TurnLin'];
kname2 = [stimnames{2} 'TurnLin'];

xyall =  btd.behaviorTriggeredDataMatrix('all', 'start', {kname1, kname2}, 0);
xyturn = btd.behaviorTriggeredDataMatrix('turn', 'start', {kname1, kname2}, 0);

turnX(:,1) = xyturn(turnvalid, 1); 
turnX(:,2) = xyturn(turnvalid, 2); 
runX(:,1) = xyall(isrun & allvalid, 1);
runX(:,2) = xyall(isrun & allvalid, 2);
allX(:,1) = xyall(allvalid, 1);
allX(:,2) = xyall(allvalid, 2);

% if(isfield(btdstruct, 'periods'))
%     
% %     n = find(strcmpi(btdstruct.periods, 'low'));
% %     theta_uv = btdstruct.theta{n};
% %     
% %     linUOp =  @(yd) cos(theta_uv).*yd{1}  + sin(theta_uv).*yd{2};
% %     linVOp = @(yd) cos(theta_uv).*yd{2} - sin(theta_uv).*yd{1};
% %     
% %     btd = btd.addOperationFields({kname1, kname2}, 'linU', linUOp);
% %     btd = btd.addOperationFields({kname1, kname2}, 'linV', linVOp );
% %     btdstruct.btd = btd;
%     
% %     uvall = btd.behaviorTriggeredDataMatrix('all', 'start', {'linU', 'linV'}, 0);
% %     uvturn = btd.behaviorTriggeredDataMatrix('turn', 'start', {'linU', 'linV'}, 0);
% 
%     turnU(:,1) = cos(theta)*turnX(:,1)+sin(theta)*turnX(:,2);
%     turnU(:,2) = -sin(theta)*turnX(:,1)+cos(theta)*turnX(:,2);
%     runU(:,1) = cos(theta)*runX(:,1)+sin(theta)*runX(:,2);
%     runU(:,2) = -sin(theta)*runX(:,1)+cos(theta)*runX(:,2);
% 
%     allU = [turnU; runU];
%     
% end

allEti = alleti(allvalid);
runEti = alleti (isrun & allvalid);
turnEti = turneti(turnvalid);
turnToff = turn_toff(turnvalid);
turnTon = turn_ton(turnvalid);
runToff = all_toff(isrun & allvalid);
runTon = all_ton(isrun & allvalid);
allToff = all_toff(allvalid);
allTon = all_ton(allvalid);
samplingInterval = median(diff(alleti(allvalid)));
qrate = btdstruct.varops.kalmanNoiseDiffusion2D;
qrate_sep = btdstruct.varops.kalmanNoiseDiffusion1D;

% for j=1:2
%     btdstruct.var_sep(j).turnX = turnX(:,j);
%     btdstruct.var_sep(j).runX = runX(:,j);
%     btdstruct.var_sep(j).allX = allX(:,j);
%     btdstruct.var_sep(j).rateEst_on = estimateRateVsTime(turnTon, turnX(:,j), runTon, runX(:,j), period, samplingInterval, qrate_sep, btdstruct.varops);
%     btdstruct.var_sep(j).rateEst_off = estimateRateVsTime(turnToff, turnX(:,j), runToff, runX(:,j), period, samplingInterval, qrate_sep, btdstruct.varops);
% end


% btdstruct.var_comb.turnX = turnX;
% btdstruct.var_comb.runX = runX;
% btdstruct.var_comb.allX = allX;
% btdstruct.var_comb.turnTon = turnTon;
% btdstruct.var_comb.turnToff = turnToff;
% btdstruct.var_comb.runTon = runTon;
% btdstruct.var_comb.runToff = runToff;
% btdstruct.var_comb.allTon = allTon;
% btdstruct.var_comb.allToff = allToff;
% btdstruct.var_comb.alleti = alleti(allvalid);

if(~isfield(btdstruct.varops, 'rateExp'))
    rateEst_off = estimateRateVsTime(turnToff, turnToff, turnX, runToff, runToff, runX, allToff, allToff, allX, period, samplingInterval, qrate, btdstruct.varops);
    for k=1:length(rateEst_off.rateExp)
        if(~isempty(rateEst_off.rateExp{k}))
            btdstruct.varops.rateExp{k}.rateParams = rateEst_off.rateExp{k}.rateParams{end};
            btdstruct.varops.rateExp{k}.theta0 = rateEst_off.rateExp{k}.theta0;
        end
    end
end

btdstruct.PPFtoff = estimateRateVsTime(turnEti, turnToff, turnX, runEti, runToff, runX, allEti, allToff, allX, period, samplingInterval, qrate, btdstruct.varops);
for k=1:length(btdstruct.PPFtoff.rateExp)
    if(~isempty(btdstruct.PPFtoff.rateExp{k}))
        btdstruct.PPFtoff.rateExp{k}.rateParams = btdstruct.varops.rateExp{k}.rateParams;
        btdstruct.PPFtoff.rateExp{k}.theta0 = btdstruct.varops.rateExp{k}.theta0;
    end
end



end

function btdstruct = ratefunGQVariance (btdstruct, opstruct, varargin)

taxis = -(0:opstruct.kernelDt:opstruct.kernelTime);

n = find([opstruct.stim.varswitch]);
gq_switch = opstruct.stim(n).gqname;
gq_conv = strcat(gq_switch, 'TurnLin');

if(strcmpi(opstruct.stim(n).ramptype, 'square'))
    gq_varh = strcat(gq_switch, '_var_high');
    gq_varl = strcat(gq_switch, '_var_low');
elseif(strcmpi(opstruct.stim(n).ramptype, 'triangle'))
    gq_varh = strcat(gq_switch, '_var_rising');
    gq_varl = strcat(gq_switch, '_var_falling');
end

gq_ton = [gq_switch, '_var_ton'];
gq_toff = [gq_switch '_var_toff'];

if (any(btdstruct.btd.findField(gq_ton) <= 0))
    btdstruct.btd = btdstruct.btd.addVarianceGQs(gq_switch, opstruct.stim(n).ramptype);
end

if (isempty(btdstruct.varops.stim(n).period))
    btdstruct.varops.stim(n).period = max(btdstruct.btd.behaviorTriggeredDataMatrix('all',[], gq_ton, 0));
end

switch_period = btdstruct.varops.stim(n).period;

btd = btdstruct.btd;

%add no-step glt

% gq_nostep = strcat(gq, '_nostep');
% gq_conv_nostep = strcat(gq_conv, '_nostep');
%
% for i=1:length(btd)
%     ind = btd(i).findField(gq);
%     old = btd(i).glt(ind);
%     new = old;
%     new.fieldname = gq_nostep;
%     new.yData(abs(new.yData)>7*std(new.yData))=0;
%     btd(i).glt(end+1) = new;
%     btd(i) = btd(i).addConvolvedFields(gq_nostep, gq_conv_nostep, btdstruct.convkernels, btdstruct.kernelDt, 'scaleToSqr', true);
% end
%
% btdstruct.btd = btd;

aa = [btd.all];
alleti = [aa.eti];
tt = [btd.turn]';
turneti = [tt.start_eti]';
isrun = [aa.isrun]';
dt = median(diff(alleti));

allton = btd.behaviorTriggeredDataMatrix('all', 'start', gq_ton, 0);
alltoff = btd.behaviorTriggeredDataMatrix('all', 'start', gq_toff, 0);
turnton = btd.behaviorTriggeredDataMatrix('turn', 'start', gq_ton, 0);
turntoff = btd.behaviorTriggeredDataMatrix('turn', 'start', gq_toff, 0);

mintime = min(alleti(allton > 0 & alltoff > 0));
turnvalid = turneti >= mintime;
allvalid = alleti >= mintime;
runvalid = allvalid & isrun;

inrun.name = 'isrun';
inrun.validop = @(x) logical(x);

ishigh = logical(btd.behaviorTriggeredDataMatrix('turn', 'start', gq_varh, 0));
islow = logical(btd.behaviorTriggeredDataMatrix('turn', 'start', gq_varl, 0));
ishigh_all = logical(btd.behaviorTriggeredDataMatrix('all', 'start', gq_varh, 0));
islow_all = logical(btd.behaviorTriggeredDataMatrix('all', 'start', gq_varl, 0));
ishigh_all = ishigh_all(runvalid);
islow_all = islow_all(runvalid);
ishigh = ishigh(turnvalid);
islow = islow(turnvalid);

% For last (tau) seconds of each cycle

tau = 5;
t1 = switch_period/2 - tau;
t2 = switch_period/2;

allton = allton(runvalid);
alltoff = alltoff(runvalid);
turnton = turnton(turnvalid);
turntoff = turntoff(turnvalid);

high_all = find(allton>=t1 & allton<=t2);
low_all = find(alltoff>=t1 & alltoff<=t2);
high = find(turnton>=t1 & turnton<=t2);
low = find(turntoff>=t1 & turntoff<=t2);


for j = find([opstruct.stim.iswn])
    

    
    gq = opstruct.stim(j).gqname;
    gq_conv = [gq 'TurnLin'];
    
    dm_all = btdstruct.btd.behaviorTriggeredDataMatrix('turn', 'start', gq, taxis);
    dm_all = dm_all(turnvalid, :);
    dm_high = dm_all(ishigh, :);
    dm2 = dm_high; dm2(~isfinite(dm2)) = 0;
    btdstruct.varRF(j).bta_high = sum(dm2)./sum(isfinite(dm_high));
    
    dm_low = dm_all(islow, :);
    dm2 = dm_low; dm2(~isfinite(dm2)) = 0;
    btdstruct.varRF(j).bta_low = sum(dm2)./sum(isfinite(dm_low));
    
    turndata = btd.behaviorTriggeredDataMatrix('turn', 'start', gq_conv, 0);
    alldata =  btd.behaviorTriggeredDataMatrix('all', 'start', gq_conv, 0);
    turndata = turndata(turnvalid);
    alldata = alldata(runvalid);
    
    nr = size(turndata, 1);
    na = size(alldata,1);
    mu = mean(turndata);
    mu_a = mean(alldata);
    s = std(turndata);
    s_a = std(alldata);
    
    nrh = nnz(ishigh);
    nrl = nnz(islow);
    nah = nnz(ishigh_all);
    nal = nnz(islow_all);
    
    muh = mean(turndata(ishigh));
    mul = mean(turndata(islow));
    mu_ah = mean(alldata(ishigh_all));
    mu_al = mean(alldata(islow_all));
    
    sh = std(turndata(ishigh));
    sl = std(turndata(islow));
    s_ah = std(alldata(ishigh_all));
    s_al = std(alldata(islow_all));
    
    ratefun = @(ydata) nr/(na*dt)*normpdf(ydata,mu,s)./normpdf(ydata, 1.0*mu_a, s_a);
    
    ratefun_high = @(ydata) nrh/(nah*dt) * normpdf(ydata, muh, sh)./normpdf(ydata, mu_ah, s_ah);
    btdstruct.btd = btdstruct.btd.addOperationFields(gq_conv, [gq '_ratePredROG_var_high'], ratefun_high);
    
    ratefun_low = @(ydata) nrl/(nal*dt) * normpdf(ydata, mul, sl)./normpdf(ydata, mu_al, s_al);
    btdstruct.btd = btdstruct.btd.addOperationFields(gq_conv, [gq '_ratePredROG_var_low'], ratefun_low);
    
    %btdstruct.btd = btdstruct.btd.addOperationFields(gq_conv_nostep, [gq_nostep '_ratePredROG'], ratefun);
    %btdstruct.btd = btdstruct.btd.addOperationFields(gq_conv_nostep, [gq_nostep '_ratePredROG_var_high'], ratefun_high);
    %btdstruct.btd = btdstruct.btd.addOperationFields(gq_conv_nostep, [gq_nostep '_ratePredROG_var_low'], ratefun_low);
    
    
    cv_high = alldata(high_all);
    cv_low = alldata(low_all);
    cvt_high = turndata(high);
    cvt_low = turndata(low);
    
    nrh_last = nnz(high);
    nrl_last = nnz(low);
    nah_last = nnz(high_all);
    nal_last = nnz(low_all);
    
    muh_last = mean(cvt_high);
    mul_last = mean(cvt_low);
    mu_ah_last = mean(cv_high);
    mu_al_last = mean(cv_low);
    
    sh_last = std(cvt_high);
    sl_last = std(cvt_low);
    s_ah_last = std(cv_high);
    s_al_last = std(cv_low);
    
    ratefun_high_last = @(ydata) nrh_last/(nah_last*dt) * normpdf(ydata, muh_last, sh_last)./normpdf(ydata, mu_ah_last, s_ah_last);
    btdstruct.btd = btdstruct.btd.addOperationFields(gq_conv, [gq '_ratePredROG_var_high_last'], ratefun_high_last);
    
    ratefun_low_last = @(ydata) nrl_last/(nal_last*dt) * normpdf(ydata, mul_last, sl_last)./normpdf(ydata, mu_al_last, s_al_last);
    btdstruct.btd = btdstruct.btd.addOperationFields(gq_conv, [gq '_ratePredROG_var_low_last'], ratefun_low_last);
    
    lx_h_last = linspace(percentile(alldata(high_all),0.01), percentile(alldata(high_all),.99), opstruct.numLxBins);
    lx_l_last = linspace(percentile(alldata(low_all),0.01), percentile(alldata(low_all),.99), opstruct.numLxBins);
    lx_h = linspace(percentile(alldata(ishigh_all),0.01), percentile(alldata(ishigh_all),.99), opstruct.numLxBins);
    lx_l = linspace(percentile(alldata(islow_all),0.01), percentile(alldata(islow_all),.99), opstruct.numLxBins);
    
    h_high = histc(alldata(ishigh_all), binEdgesFromCenters(lx_h));
    h_low = histc(alldata(islow_all), binEdgesFromCenters(lx_l));
    ht_high = histc(turndata(ishigh), binEdgesFromCenters(lx_h));
    ht_low = histc(turndata(islow), binEdgesFromCenters(lx_l));
    h_high_last = histc(alldata(high_all), binEdgesFromCenters(lx_h_last));
    h_low_last = histc(alldata(low_all), binEdgesFromCenters(lx_l_last));
    ht_high_last = histc(turndata(high), binEdgesFromCenters(lx_h_last));
    ht_low_last = histc(turndata(low), binEdgesFromCenters(lx_l_last));
    
    btdstruct.varRF(j).gqname = opstruct.stim(j).gqname;
    
    btdstruct.varRF(j).xh = lx_h;
    btdstruct.varRF(j).xh_last = lx_h_last;
    btdstruct.varRF(j).xl = lx_l;
    btdstruct.varRF(j).xl_last = lx_l_last;
    
    btdstruct.varRF(j).rog_params_high = [nrh/(nah*dt) muh, sh, mu_ah, s_ah];
    btdstruct.varRF(j).rog_params_low = [nrl/(nal*dt) mul, sl, mu_al, s_al];
    btdstruct.varRF(j).rog_params_high_last = [nrh_last/(nah_last*dt), muh_last, sh_last, mu_ah_last, s_ah_last];
    btdstruct.varRF(j).rog_params_low_last = [nrl_last/(nal_last*dt), mul_last, sl_last, mu_al_last, s_al_last];
    
    btdstruct.varRF(j).rate_rog_high =  60*ratefun_high(lx_h);
    btdstruct.varRF(j).rate_rog_low =  60*ratefun_low(lx_l);
    btdstruct.varRF(j).rate_rog_high_last =  60*ratefun_high_last(lx_h_last);
    btdstruct.varRF(j).rate_rog_low_last =  60*ratefun_low_last(lx_l_last);
    
    btdstruct.varRF(j).rate_high = 60*ht_high(1:end-1)./ h_high(1:end-1) / dt;
    btdstruct.varRF(j).rate_low = 60*ht_low(1:end-1)./ h_low(1:end-1) / dt;
    btdstruct.varRF(j).rate_high_last = 60*ht_high_last(1:end-1)./ h_high_last(1:end-1) / dt;
    btdstruct.varRF(j).rate_low_last = 60*ht_low_last(1:end-1)./ h_low_last(1:end-1) / dt;
    
end


end



