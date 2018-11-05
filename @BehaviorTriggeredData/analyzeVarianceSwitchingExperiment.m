function btdstruct = analyzeVarianceSwitchingExperiment(btdstruct, opstruct)
%btdstruct = analyzeVarianceSwitchingExperiment(btdstruct, opstruct)
%

if (nargin < 1)
    opstruct.stim(1).gqname = 'led1ValDiff';
    opstruct.stim(1).iswn = true;
    opstruct.stim(1).ramptype = 'square';
    opstruct.stim(1).switch = true;
    opstruct.stim(1).period = [];
     opstruct.stim(2).gqname = 'led2ValDiff';
     opstruct.stim(2).iswn = true;
     opstruct.stim(2).ramptype = 'square';
     opstruct.stim(2).switch = false;
     opstruct.stim(2).period = [];
    opstruct.timeRange = [-Inf Inf];
    opstruct.redogq = true;
    opstruct.kernelTime = 7;
    opstruct.kernelDt = 0.1;
    opstruct.numTimeBins = 120; %time bin size = period/numTimeBins
    opstruct.redogq = false;
    opstruct.numLxBins = 20;
    opstruct.adaptationTime = 10; %first T seconds to exclude from high/low cycles in analysis
    opstruct.Q_alpha = 0.2;
    opstruct.pdegree = 1;
    btdstruct = opstruct;
    return;
end

if (nargin < 2 || ~isfield(btdstruct, 'btd'))
    error ('first argument must be a btdstruct and second argument contains options');
end
if (~isfield(opstruct, 'pdegree'))
    opstruct.pdegree = 1;
end

if (~isfield(btdstruct, 'var') || opstruct.redogq)
    btdstruct = BehaviorTriggeredData.prepVarianceSwitchingAnalysis(btdstruct, opstruct);
else
    fn = {'stim', 'numLxBins', 'adaptationTime', 'Q_alpha', 'numTimeBins', 'timeRange'};
    for j = 1:length(fn)
        btdstruct.varops.(fn{j}) = opstruct.(fn{j});
    end
end
for j = 1:length(opstruct.stim)
    inds(j) = find(strcmpi(opstruct.stim(j).gqname, {btdstruct.var.gqname}));
end

btdstruct.var = btdstruct.var(inds);
vi = [opstruct.stim.iswn];


if (any(vi))
     btdstruct = measureRatesVsVar (btdstruct);

    %fit stretch factor vs. time
    pd = BehaviorTriggeredData.createProblemStructForRateFunFitting (btdstruct.var(vi), 'eti', opstruct.pdegree, opstruct.timeRange);
    pd.Q_alpha = opstruct.Q_alpha * eye(nnz(vi));
    pd.maxreps = 2;
    btdstruct.alphavseti = fitRateFunWithTemporalScalingND(pd);

    %turn time estimate into time in cycle estimate
    [btdstruct.alphavseti.alpha_ton, btdstruct.alphavseti.valpha_ton, btdstruct.alphavseti.tx_ton] = ...
        etiThetaToCycleTheta(btdstruct.alphavseti.alpha, btdstruct.alphavseti.valpha, btdstruct.alphavseti.tx, pd.tshift, pd.period);
    [btdstruct.alphavseti.alpha_s_ton, btdstruct.alphavseti.valpha_s_ton] = ...
        etiThetaToCycleTheta(btdstruct.alphavseti.alpha_s, btdstruct.alphavseti.valpha_s, btdstruct.alphavseti.tx, pd.tshift, pd.period);
    
    %directly calculate stretch factor vs cycle time
    
    pd = BehaviorTriggeredData.createProblemStructForRateFunFitting (btdstruct.var(vi), 'ton', 1, opstruct.timeRange);
    
    if all(size(pd.alpha_0) == size(btdstruct.alphavseti.alpha_ton))
        pd.alpha_0 = btdstruct.alphavseti.alpha_ton;
    end
    pd.params_0 = btdstruct.alphavseti.staticParams;
    pd.tparams_0 = btdstruct.alphavseti.temporalParams;
    
    pd.Q_alpha = opstruct.Q_alpha * eye(nnz(vi)) * pd.period / diff(pd.timeRange);
    [btdstruct.alphavston, scaledValues] = fitRateFunWithTemporalScalingND(pd);
    btdstruct = adaptedRatesVsVar (btdstruct, scaledValues, pd.period);
   
end
    
function btdstruct = adaptedRatesVsVar (btdstruct, scaledValues, period)
opstruct = btdstruct.varops;
tton = scaledValues.turnT;
teti = scaledValues.turnEti;
rton = scaledValues.runT;
reti = scaledValues.runEti;

tval = scaledValues.tval;
rval = scaledValues.rval;
deltaT = scaledValues.deltaT;

tvalid = teti >= min(opstruct.timeRange) & teti < max(opstruct.timeRange);
thigh = tton > opstruct.adaptationTime & tton < period/2 & tvalid;
tlow = tton > period/2 + opstruct.adaptationTime & tton < period & tvalid;
tup = tton > 0 & tton < opstruct.adaptationTime & tvalid;
tdown = tton > period/2  & tton < period/2 + opstruct.adaptationTime & tvalid;


rvalid = reti >= min(opstruct.timeRange) & reti < max(opstruct.timeRange);
rhigh = rton > opstruct.adaptationTime & rton < period/2 & rvalid;
rlow = rton > period/2 + opstruct.adaptationTime & rton < period & rvalid;
rup = rton > 0 & rton < opstruct.adaptationTime & rvalid;
rdown = rton > period/2  & rton < period/2 + opstruct.adaptationTime & rvalid;

 

for j = 1:size(rval, 2)
    lx = linspace(-3*std(rval(rvalid,j)), 3*std(rval(rvalid,j)), opstruct.numLxBins);
    btdstruct.scaledRFHigh(j) = turnRates(tval(thigh, j), rval(rhigh, j), deltaT, lx);
    btdstruct.scaledRFLow(j) = turnRates(tval(tlow, j), rval(rlow, j), deltaT, lx);
    btdstruct.scaledRFUp(j) = turnRates(tval(tup, j), rval(rup, j), deltaT, lx);
    btdstruct.scaledRFDown(j) = turnRates(tval(tdown, j), rval(rdown, j), deltaT, lx);
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


    

function btdstruct = measureRatesVsVar (btdstruct)

opstruct = btdstruct.varops;
vi = find([opstruct.stim.iswn]);
if (isempty(vi))
    return;
end

op = optimoptions('fminunc');
op.Display = 'off';
op.Algorithm = 'quasi-newton';

for j = vi
    deltaT = median(diff(btdstruct.var(j).fullensemble.eti));
  
    ti = btdstruct.var(j).turn.ton > opstruct.adaptationTime & btdstruct.var(j).turn.ton < btdstruct.var(j).turn.toff & ...
        btdstruct.var(j).turn.eti >= min(opstruct.timeRange) & btdstruct.var(j).turn.eti < max(opstruct.timeRange);
    ri = btdstruct.var(j).noturn.ton > opstruct.adaptationTime & btdstruct.var(j).noturn.ton < btdstruct.var(j).noturn.toff & ...
        btdstruct.var(j).noturn.eti >= min(opstruct.timeRange) & btdstruct.var(j).noturn.eti < max(opstruct.timeRange);
    
    xt = btdstruct.var(j).turn.x_conv(ti);
    xr = btdstruct.var(j).noturn.x_conv(ri);
    
    
    lx = linspace(-3*std(xr), 3*std(xr), opstruct.numLxBins);
    
    btdstruct.rfHighVar(j) = turnRates (xt, xr, deltaT, lx);
    
    
    ti = btdstruct.var(j).turn.toff > opstruct.adaptationTime & btdstruct.var(j).turn.toff < btdstruct.var(j).turn.ton & ...
        btdstruct.var(j).turn.eti >= min(opstruct.timeRange) & btdstruct.var(j).turn.eti < max(opstruct.timeRange);
    ri = btdstruct.var(j).noturn.toff > opstruct.adaptationTime & btdstruct.var(j).noturn.toff < btdstruct.var(j).noturn.ton & ...
        btdstruct.var(j).noturn.eti >= min(opstruct.timeRange) & btdstruct.var(j).noturn.eti < max(opstruct.timeRange);
    
    xt = btdstruct.var(j).turn.x_conv(ti);
    xr = btdstruct.var(j).noturn.x_conv(ri);
    
    lx = linspace(-3*std(xr), 3*std(xr), opstruct.numLxBins);
    
    
    btdstruct.rfLowVar(j)  = turnRates (xt, xr, deltaT, lx);
   
end
    


