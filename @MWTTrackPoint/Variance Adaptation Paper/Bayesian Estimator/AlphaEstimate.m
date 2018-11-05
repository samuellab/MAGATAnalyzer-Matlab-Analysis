
function out = AlphaEstimate(btdstruct, opstruct, timeType, isconvolved, Ddt, deltaT, taxis, period, tshift)
% function [btdstruct, bayesalpha] = AlphaEstimate(btdstruct, opstruct, RescaleStruct, isconvolved, Ddt, deltaT, taxis, period, tshift)
% Finds the stimulus at deltaT intervals, and calculates the Bayes-Optimal estimates of its variance


if(~isfield(btdstruct, 'btd'))
    disp('cannot calculate Bayes-Optimal estimate without btds');
    return
end

timeField = 'eti';
nexps = length(btdstruct.btd);
switchStim = find([opstruct.stim.switch]);
gq = opstruct.stim(switchStim).gqname;
    
gq_ton = [gq '_var_ton'];
gq_toff = [gq '_var_toff'];

if(strcmpi(gq, 'led1ValDiff'))
    ind = findField(btdstruct.btd(1), 'led1Val');
    gq_conv = [gq 'TurnLin'];
elseif(strcmpi(gq, 'led2ValDiff'))
    ind = findField(btdstruct.btd(1), 'led2Val');
    gq_conv = [gq 'TurnLin'];
elseif(strcmpi(gq, 'led2m1ValDiff'))
    ind = findField(btdstruct.btd(1), 'led2m1Val');
    gq_conv = 'linU';
end

if(~isfield(btdstruct, 'var') || ind<=0 || any(btdstruct.btd.findField(gq_conv) <= 0))
    btdstruct = BehaviorTriggeredData.prepVarianceSwitchingAnalysis_Gepner(btdstruct, opstruct);
end
if(nargin<7)
    pd = BehaviorTriggeredData.createProblemStructForRateFunFitting (btdstruct.var, timeField, opstruct.pdegree, opstruct.timeRange, [], deltaT);
    taxis = pd.tx;
    period = pd.period;
    tshift = pd.tshift;
end



if(strcmpi(gq, 'led1ValDiff'))
    ind = findField(btdstruct.btd(1), 'led1Val');
    gq_conv = [gq 'TurnLin'];
elseif(strcmpi(gq, 'led2ValDiff'))
    ind = findField(btdstruct.btd(1), 'led2Val');
    gq_conv = [gq 'TurnLin'];
elseif(strcmpi(gq, 'led2m1ValDiff'))
    ind = findField(btdstruct.btd(1), 'led2m1Val');
    gq_conv = 'linU';
end

ind_diff = findField(btdstruct.btd(1), gq);
ind_conv = findField(btdstruct.btd(1), gq_conv);
ind_on = findField(btdstruct.btd(1), gq_ton);

btd = btdstruct.btd;

mS_all_eti = [];
vS_all_eti = [];
mS_all_ton = [];
vS_all_ton = [];

tx = mod(taxis+tshift, period);
tcycle = 0:deltaT:period;

for i=1:nexps
    
    glt = btd(i).glt(ind);
    glt_conv = btd(i).glt(ind_conv);
    glt_on = btd(i).glt(ind_on);
    glt_diff = btd(i).glt(ind_diff);
    
    xval = interp1(glt.xData, glt.yData, taxis); % stimulus
    xd = diff(xval);
    xdiff(i,:) = [0 xd]; % stimulus derivative
    
    xdiff2(i,:) = interp1(glt_diff.xData, glt_diff.yData/120, taxis);
    
    xo(i,:) = interp1(glt_conv.xData, glt_conv.yData, taxis); % convolved stimulus
    

    
    [~,~,~,sconv] = meanyvsx(tx, xo(i,:), tcycle);
    [~,~,~,s] = meanyvsx(tx, xdiff(i,:), tcycle);
    
    xdiff(i,:) = xdiff(i,:)./mean(s(round(end/2)+1:end)); % normalize the amplitude of the white noise stimulus derivative
    xo(i,:) = xo(i,:)./mean(sconv(round(end/2)+1:end));
    
    [~,~,~,son1] = meanyvsx(tx, xdiff(i,:), tcycle);
%     [~,~,~,son2] = meanyvsx(tx, xdiff2(i,:), tcycle);
    [~,~,~,son_conv] = meanyvsx(tx, xo(i,:), tcycle);
    
    
    if(isconvolved)
        [mS_eti(i,:), vS_eti(i,:)] = BayesianVarEstimate(1, Ddt, xo(i,:));
    else
        [mS_eti(i,:), vS_eti(i,:)] = BayesianVarEstimate(1, Ddt, xdiff(i,:));
    end
    [mS_ton(i,:), vS_ton(i,:), tx_ton] = meanEstinCycle(mS_eti(i,:), vS_eti(i,:), taxis, tshift, period, opstruct.timeRange, timeType);
    
    mS_all_eti = [mS_all_eti mS_eti(i,:)];
    vS_all_eti = [vS_all_eti vS_eti(i,:)];
    mS_all_ton = [mS_all_ton mS_ton(i,:)];
    vS_all_ton = [vS_all_ton vS_ton(i,:)];

end

[msigma, vsigma, tt] = meanEstinCycle(mS_all_eti, vS_all_eti, taxis, 0, max(taxis), opstruct.timeRange, 'ton');
[msigma_ton, vsigma_ton, tton] = meanEstinCycle(mS_all_ton, vS_all_ton, tx_ton, tshift, period, opstruct.timeRange, 'ton');

out.tx = tt;
out.sigma = msigma;
out.vsigma = vsigma;
out.tx_ton = tton;
out.sigma_ton = msigma_ton;
out.vsigma_ton = vsigma_ton;

end