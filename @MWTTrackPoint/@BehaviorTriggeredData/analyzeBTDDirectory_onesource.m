function btdstruct = analyzeBTDDirectory_onesource (btdstruct, gqname, varargin)
%function btdstruct = analyzeBTDDirectory_onesource (btdstruct, gqname, varargin)


taxis = -30:0.1:10;
kernelTime = 7; %s
kernelDt = 0.1;
turnSizeRange = [-15 -3];
headSweepAcceptRange = [0 1.25];
compensatorName = '';
valName = 'led1Val';
%deltaTPower = 1; %derivatives
varargin = assignApplicable(varargin);



if (~isstruct(btdstruct))
    if (ischar(btdstruct) || (iscell(btdstruct) && all(cellfun(@ischar, btdstruct))))
        btdstruct = BehaviorTriggeredData.loadBTDDirectory(btdstruct);
    else
        error ('first argument must be btdstruct or basedir');
    end
end
btd = btdstruct.btd;
if (isempty(btd))
    return;
end

if (~existsAndDefault('gqname', 'led1ValDiff')) 
    if (isfield(btdstruct, 'gqname'))
        gqname = btdstruct.gqname;
    end
end
%existsAndDefault('compensatorName', '');

btdstruct.kernelTime = kernelTime;
btdstruct.kernelDt = kernelDt;
btdstruct.taxis = taxis;
btdstruct.turnSizeRange = turnSizeRange;
btdstruct.headSweepAcceptRange = headSweepAcceptRange;
btdstruct.gqname = gqname;
hs = [btd.hs];
hsd = [hs.end_eti] - [hs.start_eti];
hsa = logical([hs.accepted]);
btdstruct.hsdur.mean = mean(hsd);
btdstruct.hsdur.median = median(hsd);
btdstruct.hsdur.sigma = std(hsd);

btdstruct.hsdur.mean_acc = mean(hsd(hsa));
btdstruct.hsdur.median_acc = median(hsd(hsa));
btdstruct.hsdur.sigma_acc = std(hsd(hsa));

btdstruct.hsdur.mean_rej = mean(hsd(~hsa));
btdstruct.hsdur.median_rej = median(hsd(~hsa));
btdstruct.hsdur.sigma_rej = std(hsd(~hsa));


for j = 1:length(btd)
    btd(j).all.absdeltatheta = abs(btd(j).all.deltatheta);
    tn = btd(j).all.tnum;
    btd(j).all.absdeltatheta = abs(btd(j).all.deltatheta);
    btd(j).all.dabsdeltatheta = zeros(size(btd(j).all.absdeltatheta));
    for k = 1:max(tn)
        btd(j).all.dabsdeltatheta(tn == k) = [0 diff(btd(j).all.absdeltatheta(tn == k))];
    end
    btd(j).all.absspineTheta = abs(btd(j).all.spineTheta);
    btd(j).all.dabsspineTheta = zeros(size(btd(j).all.absspineTheta));
    for k = 1:max(tn)
        btd(j).all.dabsspineTheta(tn == k) = [0 diff(btd(j).all.absspineTheta(tn == k))];
    end
    btd(j).all.abscurv = abs(btd(j).all.curv);
    pct = percentile(btd(j).all.abscurv, 0.95);
    btd(j).all.abscurv = min(btd(j).all.abscurv, pct);
end

fields = {'speed', 'absdeltatheta','dabsdeltatheta', 'abscurv', 'absspineTheta','dabsspineTheta'};
for j = 1:length(fields)
    [k0,k1] = btd.wienerKernels(gqname, fields{j}, -taxis);
    btdstruct.wiener.(fields{j}).k0 = k0;
    btdstruct.wiener.(fields{j}).k1 = k1;
end


fields = {'turn', 'pause', 'acchs', 'rejhs'};% 'acchs', 'rejhs','acchs', 'rejhs'};
pos =   {'start', 'start', 'start', 'start'};% 'atmax', 'atmax',  'end',   'end'};

for j = 1:length(fields)
    [bta, ~, eb] = btd.behaviorTriggeredAverage(fields{j}, pos{j}, gqname, taxis);
    btdstruct.([fields{j} '_' pos{j}]) = bta;
    btdstruct.([fields{j} '_' pos{j} '_eb']) = eb;
    if (~isempty(compensatorName))
        btdstruct.([fields{j} '_' pos{j} '_compensator']) = btd.behaviorTriggeredAverageCompensator(fields{j}, pos{j}, gqname, compensatorName, taxis);
    end
    try
        [bta, ~, eb] = btd.behaviorTriggeredAverage(fields{j}, pos{j}, valName, taxis);
        btdstruct.val.([fields{j} '_' pos{j}]) = bta;
        btdstruct.val.([fields{j} '_' pos{j} '_eb']) = eb;
    catch me %#ok<NASGU>
    end
    
end

fields = {'turn'};% 'acchs', 'rejhs','acchs', 'rejhs'};
pos =   {'start'};% 'atmax', 'atmax',  'end',   'end'};

for k = 1:length(btd)
    
    for j = 1:length(fields)
        btdstruct.exp(k).([fields{j} '_' pos{j}]) = btd(k).behaviorTriggeredAverage(fields{j}, pos{j}, gqname, taxis)';
    end
end


kname = [gqname 'TurnLin'];
btdstruct.linname = kname;

inrun.name = 'isrun';
inrun.validop = @(x) logical(x);

[btdstruct.convkernels, btd] = btd.createBTAKernel('turn', 'start', gqname, kernelTime, kernelDt, 'newFieldName',kname, 'abbott', true);
s = std(btd.behaviorTriggeredDataMatrix ('all', 'start', kname, 0));
% btdstruct.convkernels = btdstruct.convkernels / s;
btd = btd.addConvolvedFields(gqname, kname, btdstruct.convkernels, kernelDt, 'scaleToSqr', true);

[btdstruct.tr_vs_conv, btdstruct.tr_vs_conv_eb, btdstruct.lx] = rateVsField(btd, 'turn', 'start', kname);
for k = 1:length(btd)
    [tr_vs_conv, tr_vs_conv_eb, lx] = rateVsField(btd(k), 'turn', 'start', kname);
    btdstruct.exp(k).tr_vs_conv = tr_vs_conv';
    btdstruct.exp(k).tr_vs_conv_eb = tr_vs_conv_eb';
    btdstruct.exp(k).lx = lx';
end



firstten.name = 'start_eti';
firstten.validop = @(x) x < 600;
secondten = firstten;
secondten.validop = @(x) (x >= 600 & x < 1200);

fields = {'turn', 'acchs', 'rejhs'};% 'acchs', 'rejhs','acchs', 'rejhs'};
pos =   {'start', 'start', 'start'};% 'atmax', 'atmax',  'end',   'end'};

for j = 1:length(fields)
    [bta, ~, eb] = btd.behaviorTriggeredAverage(fields{j}, pos{j}, gqname, taxis, 'conditions', firstten);
    btdstruct.firstten.([fields{j} '_' pos{j}]) = bta;
    btdstruct.firstten.([fields{j} '_' pos{j} '_eb']) = eb;
    [bta, ~, eb] = btd.behaviorTriggeredAverage(fields{j}, pos{j}, gqname, taxis, 'conditions', secondten);
    btdstruct.secondten.([fields{j} '_' pos{j}]) = bta;
    btdstruct.secondten.([fields{j} '_' pos{j} '_eb']) = eb;
   
end


fieldlist = {kname};
ratenames = {[gqname '_ratePredROG']};
turndata = btd.behaviorTriggeredDataMatrix('turn', 'start', fieldlist, 0);
alldata =  btd.behaviorTriggeredDataMatrix('all', 'start', fieldlist, 0, 'conditions', inrun);

tt = [btd.turn];
turneti = [tt.start_eti];
aa = [btd.all];
alleti = [aa.eti]; alleti = alleti(logical([aa.isrun]));

dt = median(diff(btd(1).all.eti)); %assuming all movies were captured at same frame rate!


xx = binEdgesFromCenters(btdstruct.lx);
for j = 1:length(fieldlist)
    yy = turndata(turneti < 600, j);
    aa = alldata(alleti < 600, j);
    h1 = histc(yy,xx);
    h2 = histc(aa,xx);
    rr = 1/dt*h1(1:(end-1))./h2(1:(end-1));
    reb =  1/dt*sqrt(h1(1:(end-1)))./h2(1:(end-1));
    btdstruct.firstten.tr_vs_conv(j,:) = rr';
    btdstruct.firstten.tr_vs_conv_eb(j,:) = reb';
    
end
for j = 1:length(fieldlist)
    yy = turndata(turneti >= 600 & turneti < 1200, j);
    aa = alldata(alleti >= 600 & alleti < 1200, j);
    h1 = histc(yy,xx);
    h2 = histc(aa,xx);
    rr = 1/dt*h1(1:(end-1))./h2(1:(end-1));
    reb =  1/dt*sqrt(h1(1:(end-1)))./h2(1:(end-1));
    btdstruct.secondten.tr_vs_conv(j,:) = rr';
    btdstruct.secondten.tr_vs_conv_eb(j,:) = reb';
end


btdstruct.firstten.nr = nnz(turneti < 600);
btdstruct.firstten.na = nnz(alleti<600);
btdstruct.firstten.mu = mean(turndata(turneti < 600,:));
btdstruct.firstten.mu_a = mean(alldata(alleti < 600,:));
btdstruct.firstten.s = std(turndata(turneti < 600,:));
btdstruct.firstten.s_a = std(alldata(alleti < 600,:));

btdstruct.secondten.nr = nnz(turneti >= 600 & turneti < 1200);
btdstruct.secondten.na = nnz(alleti >= 600 & alleti < 1200);
btdstruct.secondten.mu = mean(turndata(turneti >= 600 & turneti < 1200,:));
btdstruct.secondten.mu_a = mean(alldata(alleti >= 600 & alleti < 1200,:));
btdstruct.secondten.s = std(turndata(turneti >= 600 & turneti < 1200,:));
btdstruct.secondten.s_a = std(alldata(alleti >= 600 & alleti < 1200,:));


nr = size(turndata, 1);
na = size(alldata,1);
mu = mean(turndata);
mu_a = mean(alldata);
s = std(turndata);
s_a = std(alldata);

for j = 1:size(turndata, 2)
    
    ratefun = @(ydata) nr/(na*dt)*normpdf(ydata,mu(j),s(j))./normpdf(ydata, 1.0*mu_a(j), s_a(j));
    btdstruct.ratefun{j} = ratefun;
    btd = btd.addOperationFields(fieldlist{j}, ratenames{j}, ratefun);
    
    %{
    ratefun0 = @(x, xdata) x(1) * normpdf(xdata, x(2), x(3)) ./ normpdf(xdata, 0, x(4)); 
    initguess0 = [nr/(na*dt) mu(j) s(j) s_a(j)];
    fitfun0 = @(x) -sum(log(max(ratefun0(x,turndata(:,j)),1E-100))) + sum(max(0,ratefun0(x,alldata(:,j))));
    
    op = optimset('fmincon');
    op.Algorithm = 'active-set';
    %op.LargeScale = 'off';
    problem.objective = fitfun0;
    problem.solver = 'fmincon';
    problem.x0 = initguess0;
    problem.lb = [0 -Inf 0  0];
    problem.ub = [Inf Inf Inf  Inf];
    problem.options = op;
    fitParams0 = fmincon(problem);
    
    ratefun = @(ydata) ratefun0(fitParams0, ydata);
    btd = btd.addOperationFields(fieldlist{j}, [ratenames{j} '_fit'], ratefun);
    %}

    
    
    ratefun = @(ydata) btdstruct.firstten.nr/(btdstruct.firstten.na*dt)*normpdf(ydata,btdstruct.firstten.mu(j),btdstruct.firstten.s(j))...
        ./normpdf(ydata, 1.0*btdstruct.firstten.mu_a(j), btdstruct.firstten.s_a(j));
    btdstruct.firstten.ratefun{j} = ratefun;
    btd = btd.addOperationFields(fieldlist{j}, [ratenames{j} '_firstten'], ratefun);
        
    ratefun = @(ydata) btdstruct.secondten.nr/(btdstruct.secondten.na*dt)*normpdf(ydata,btdstruct.secondten.mu(j),btdstruct.secondten.s(j))...
        ./normpdf(ydata, 1.0*btdstruct.secondten.mu_a(j), btdstruct.secondten.s_a(j));
    btdstruct.secondten.ratefun{j} = ratefun;
    btd = btd.addOperationFields(fieldlist{j}, [ratenames{j} '_secondten'], ratefun);
end

btdstruct.rate_rog =  btdstruct.ratefun{1}(btdstruct.lx);
btdstruct.firstten.rate_rog =  btdstruct.firstten.ratefun{1}(btdstruct.lx);
btdstruct.secondten.rate_rog =  btdstruct.secondten.ratefun{1}(btdstruct.lx);


ratenames_firstten = cellfun(@(s) [s '_firstten'], ratenames, 'UniformOutput', false);
ratenames_secondten = cellfun(@(s) [s '_secondten'], ratenames, 'UniformOutput', false);
btdstruct.predratefields = [ratenames ratenames_firstten ratenames_secondten];
turn_predrate = btd.behaviorTriggeredDataMatrix('turn', 'start', btdstruct.predratefields, 0);
all_predrate =  btd.behaviorTriggeredDataMatrix('all', 'start', btdstruct.predratefields, 0, 'conditions', inrun);

btdstruct.trx = linspace(min(btdstruct.ratefun{1}(btdstruct.lx)), max(btdstruct.ratefun{1}(btdstruct.lx)), length(btdstruct.lx));

trxx = binEdgesFromCenters(btdstruct.trx);
btdstruct.trpred_rsq = [];
for j = 1:length(btdstruct.predratefields)
    yy = turn_predrate(:, j);
    aa = all_predrate(:, j);
    h1 = histc(yy,trxx);
    h2 = histc(aa,trxx);
    rr = 1/dt*h1(1:(end-1))./h2(1:(end-1));
    reb =  1/dt*sqrt(h1(1:(end-1)))./h2(1:(end-1));
    btdstruct.tr_vs_pred(j,:) = rr';
    btdstruct.tr_vs_pred_eb(j,:) = reb'; 
    btdstruct.nturns_vs_pred(j,:) = h1(1:(end-1))';
    mtr = sum(btdstruct.nturns_vs_pred(j,:).*btdstruct.tr_vs_pred(j,:))./sum(btdstruct.nturns_vs_pred(j,:));
    sstot =  sum(btdstruct.nturns_vs_pred(j,:).*(btdstruct.tr_vs_pred(j,:)-mtr).^2)./sum(btdstruct.nturns_vs_pred(j,:));
    ssres =  sum(btdstruct.nturns_vs_pred(j,:).*(btdstruct.tr_vs_pred(j,:)-btdstruct.trx).^2)./sum(btdstruct.nturns_vs_pred(j,:));
    btdstruct.trpred_rsq(j) = 1 - ssres./sstot;
    
end
btdstruct.firstten.trpred_rsq = [];
for j = 1:length(btdstruct.predratefields)
    yy = turn_predrate(turneti < 600, j);
    aa = all_predrate(alleti < 600, j);
    h1 = histc(yy,trxx);
    h2 = histc(aa,trxx);
    rr = 1/dt*h1(1:(end-1))./h2(1:(end-1));
    reb =  1/dt*sqrt(h1(1:(end-1)))./h2(1:(end-1));
    btdstruct.firstten.tr_vs_pred(j,:) = rr';
    btdstruct.firstten.tr_vs_pred_eb(j,:) = reb';
    btdstruct.firstten.nturns_vs_pred(j,:) = h1(1:(end-1))';
    mtr = sum(btdstruct.firstten.nturns_vs_pred(j,:).*btdstruct.firstten.tr_vs_pred(j,:))./sum(btdstruct.firstten.nturns_vs_pred(j,:));
    sstot =  sum(btdstruct.firstten.nturns_vs_pred(j,:).*(btdstruct.firstten.tr_vs_pred(j,:)-mtr).^2)./sum(btdstruct.firstten.nturns_vs_pred(j,:));
    ssres =  sum(btdstruct.firstten.nturns_vs_pred(j,:).*(btdstruct.firstten.tr_vs_pred(j,:)-btdstruct.trx).^2)./sum(btdstruct.firstten.nturns_vs_pred(j,:));
    btdstruct.firstten.trpred_rsq(j) = 1 - ssres./sstot;
    
end
btdstruct.secondten.trpred_rsq = [];
for j = 1:length(btdstruct.predratefields)
    yy = turn_predrate(turneti >= 600 & turneti < 1200, j);
    aa = all_predrate(alleti >= 600 & alleti < 1200, j);
    h1 = histc(yy,trxx);
    h2 = histc(aa,trxx);
    rr = 1/dt*h1(1:(end-1))./h2(1:(end-1));
    reb =  1/dt*sqrt(h1(1:(end-1)))./h2(1:(end-1));
    btdstruct.secondten.tr_vs_pred(j,:) = rr';
    btdstruct.secondten.tr_vs_pred_eb(j,:) = reb';
    btdstruct.secondten.nturns_vs_pred(j,:) = h1(1:(end-1))';
    mtr = sum(btdstruct.secondten.nturns_vs_pred(j,:).*btdstruct.secondten.tr_vs_pred(j,:))./sum(btdstruct.secondten.nturns_vs_pred(j,:));
    sstot =  sum(btdstruct.secondten.nturns_vs_pred(j,:).*(btdstruct.secondten.tr_vs_pred(j,:)-mtr).^2)./sum(btdstruct.secondten.nturns_vs_pred(j,:));
    ssres =  sum(btdstruct.secondten.nturns_vs_pred(j,:).*(btdstruct.secondten.tr_vs_pred(j,:)-btdstruct.trx).^2)./sum(btdstruct.secondten.nturns_vs_pred(j,:));
    btdstruct.secondten.trpred_rsq(j) = 1 - ssres./sstot;
    
end




cv = btd.behaviorTriggeredDataMatrix('all', 'start', btdstruct.linname, 0, 'conditions', inrun);
cvt = btd.behaviorTriggeredDataMatrix('turn', 'start', btdstruct.linname, 0);
cvt_first = btd.behaviorTriggeredDataMatrix('turn', 'start', btdstruct.linname, 0, 'conditions', firstten);
cvt_second = btd.behaviorTriggeredDataMatrix('turn', 'start', btdstruct.linname, 0, 'conditions', secondten);

btdstruct.hist_x = linspace(percentile(cv,0.01), percentile(cv,.99), 20);

btdstruct.turn_linval = cvt;
btdstruct.firstten.turn_linval = cvt_first;
btdstruct.secondten.turn_linval = cvt_second;
hx = btdstruct.hist_x;
h1 = histc(cv, binEdgesFromCenters(hx));
btdstruct.hist_all_val = h1(1:end-1);
h2 = histc(cvt, binEdgesFromCenters(hx));
btdstruct.hist_turn_val = h2(1:end-1);

btdstruct.rate = 60*btdstruct.hist_turn_val./btdstruct.hist_all_val/dt;

btdstruct.linval_mean_all = mean(cv);
btdstruct.linval_mean_turn = mean(cvt);
btdstruct.linval_var_all = var(cv);
btdstruct.linval_var_turn = var(cvt);



m1 = mean(cv);
s1 = std(cv);

cvt = (cvt-m1)/s1;
btdstruct.dkl_est_gauss = 0.5*(var(cvt) - log(var(cvt)) + mean(cvt).^2 - 1);
%btdstruct.dkl_est_gauss = 0.5*(var(cvt)/var(cv) - log(var(cvt)/var(cv)) + mean(cvt.^2)/var(cv)-1);


p = polyfit(btdstruct.lx, dt*btdstruct.tr_vs_conv,4);
btdstruct.turnratepoly = fitRateFun(btd, 'turn', 'start', kname, @polyval, p);

btd = btd.addOperationFields (kname, [gqname '_ratePredPoly'], @(x) 1/dt*polyval(btdstruct.turnratepoly, x));
firsths.name = 'hsnum';
firsths.validop = @(x) x == 1;

fields = {'acchs', 'rejhs'};%, 'acchs', 'rejhs','acchs', 'rejhs'};
pos =   {'start', 'start'};%, 'atmax', 'atmax',  'end',   'end'};

for j = 1:length(fields)
    [bta, ~] = btd.behaviorTriggeredAverage(fields{j}, pos{j}, gqname, taxis, 'conditions', firsths);
     btdstruct.(['first_' fields{j} '_' pos{j}]) = bta;
end


fields = {'acchs', 'rejhs'};%, 'acchs', 'rejhs','acchs', 'rejhs'};
pos =   {'start', 'start'};%, 'atmax', 'atmax',  'end',   'end'};

for j = 1:length(fields)
    [bta, ~] = btd.behaviorTriggeredAverage(fields{j}, pos{j}, gqname, taxis, 'conditions', [firstten firsths]);
     btdstruct.firstten.(['first_' fields{j} '_' pos{j}]) = bta;
     [bta, ~] = btd.behaviorTriggeredAverage(fields{j}, pos{j}, gqname, taxis, 'conditions', [secondten firsths]);
     btdstruct.secondten.(['first_' fields{j} '_' pos{j}]) = bta;
end

t = [btd.turn]; 
btdstruct.meanSquaredTurn = mean([t.dtheta].^2);
ms = btdstruct.meanSquaredTurn;
bigturn.name = 'dtheta';
bigturn.validop = @(x) x.^2 > ms;
smallturn.name = 'dtheta';
smallturn.validop = @(x) x.^2 < ms;
onehs.name = 'numhs';
onehs.validop = @(x) x == 1;
multihs = onehs;
multihs.validop = @(x) x > 1;
if (~isempty(compensatorName))
    btdstruct.bigturn_start_all_compensator = btd.behaviorTriggeredAverageCompensator('turn', 'start', gqname, compensatorName, taxis, 'conditions', bigturn);
    btdstruct.smallturn_start_all_compensator = btd.behaviorTriggeredAverageCompensator('turn', 'start', gqname, compensatorName, taxis, 'conditions', smallturn);
    btdstruct.smallturn_start_onehs_compensator = btd.behaviorTriggeredAverageCompensator('turn', 'start', gqname, compensatorName,taxis, 'conditions', [smallturn onehs]);
    btdstruct.bigturn_start_onehs_compensator = btd.behaviorTriggeredAverageCompensator('turn', 'start', gqname, compensatorName, taxis, 'conditions', [bigturn onehs]);
end

btdstruct.bigturn_start_onehs = btd.behaviorTriggeredAverage('turn', 'start', gqname, taxis, 'conditions', [bigturn onehs]);
btdstruct.smallturn_start_onehs = btd.behaviorTriggeredAverage('turn', 'start', gqname, taxis, 'conditions', [smallturn onehs]);
btdstruct.bigturn_start_multis = btd.behaviorTriggeredAverage('turn', 'start', gqname, taxis, 'conditions', [bigturn multihs]);
btdstruct.smallturn_start_multihs = btd.behaviorTriggeredAverage('turn', 'start', gqname, taxis, 'conditions', [smallturn multihs]);

btdstruct.bigturn_start_all = btd.behaviorTriggeredAverage('turn', 'start', gqname, taxis, 'conditions', bigturn);
btdstruct.smallturn_start_all = btd.behaviorTriggeredAverage('turn', 'start', gqname, taxis, 'conditions', smallturn);


btdstruct.weightedTurnStart_multihs = btd.behaviorTriggeredWeightedAverage('turn', 'start', 'dtheta',gqname, taxis, 'weightedOp', @abs, 'conditions', multihs);
btdstruct.weightedTurnStart_onehs = btd.behaviorTriggeredWeightedAverage('turn', 'start', 'dtheta', gqname, taxis, 'weightedOp', @abs, 'conditions', onehs);
btdstruct.weightedTurnStart_all = btd.behaviorTriggeredWeightedAverage('turn', 'start', 'dtheta', gqname, taxis, 'weightedOp', @abs);

try
    btdstruct.val.bigturn_start_onehs = btd.behaviorTriggeredAverage('turn', 'start', valName, taxis, 'conditions', [bigturn onehs]);
    btdstruct.val.smallturn_start_onehs = btd.behaviorTriggeredAverage('turn', 'start', valName, taxis, 'conditions', [smallturn onehs]);
    btdstruct.val.bigturn_start_multis = btd.behaviorTriggeredAverage('turn', 'start', valName, taxis, 'conditions', [bigturn multihs]);
    btdstruct.val.smallturn_start_multihs = btd.behaviorTriggeredAverage('turn', 'start', valName, taxis, 'conditions', [smallturn multihs]);

    btdstruct.val.bigturn_start_all = btd.behaviorTriggeredAverage('turn', 'start', valName, taxis, 'conditions', bigturn);
    btdstruct.val.smallturn_start_all = btd.behaviorTriggeredAverage('turn', 'start', valName, taxis, 'conditions', smallturn);
    
    btdstruct.val.weightedTurnStart_multihs = btd.behaviorTriggeredWeightedAverage('turn', 'start', 'dtheta',valName, taxis, 'weightedOp', @abs, 'conditions', multihs);
    btdstruct.val.weightedTurnStart_onehs = btd.behaviorTriggeredWeightedAverage('turn', 'start', 'dtheta', valName, taxis, 'weightedOp', @abs, 'conditions', onehs);
    btdstruct.val.weightedTurnStart_all = btd.behaviorTriggeredWeightedAverage('turn', 'start', 'dtheta', valName, taxis, 'weightedOp', @abs);

catch me %#ok<NASGU>
end
num_tspts = 200;

dm = btd.behaviorTriggeredDataMatrix('turn', 'start', gqname, linspace(turnSizeRange(1), turnSizeRange(2), num_tspts));
dm2 = dm; dm2(~isfinite(dm2)) = 0;
turnmc = mean(dm2, 2);

t = [btd.turn]; 
isbigturn = ([t.dtheta].^2 > btdstruct.meanSquaredTurn);
issmallturn = ~isbigturn;
hasonehs = [t.numhs] == 1;

btdstruct.meanChangeBeforeTurn = mean(turnmc, 1);
btdstruct.meanChangeBeforeTurn_eb = std(turnmc, 1)./sqrt(length(turnmc));

btdstruct.meanChangeBeforeTurn_onehs = sum(turnmc(hasonehs,:), 1)./nnz(hasonehs);
btdstruct.meanChangeBeforeTurn_onehs_eb = std(turnmc(hasonehs,:), 1)./sqrt(nnz(hasonehs));

btdstruct.meanChangeBeforeBigTurn = sum(turnmc(isbigturn,:), 1)./nnz(isbigturn);
btdstruct.meanChangeBeforeBigTurn_eb = std(turnmc(isbigturn,:), 1)./sqrt(nnz(isbigturn));

btdstruct.meanChangeBeforeSmallTurn = sum(turnmc(issmallturn,:), 1)./nnz(issmallturn);
btdstruct.meanChangeBeforeSmallTurn_eb = std(turnmc(issmallturn,:), 1)./sqrt(nnz(issmallturn));

btdstruct.meanChangeBeforeBigTurn_onehs = sum(turnmc(isbigturn & hasonehs,:), 1)./nnz(isbigturn & hasonehs);
btdstruct.meanChangeBeforeBigTurn_onehs_eb = std(turnmc(isbigturn & hasonehs,:), 1)./sqrt(nnz(isbigturn & hasonehs));

btdstruct.meanChangeBeforeSmallTurn_onehs = sum(turnmc(issmallturn & hasonehs,:), 1)./nnz(issmallturn & hasonehs);
btdstruct.meanChangeBeforeSmallTurn_onehs_eb = std(turnmc(issmallturn & hasonehs,:), 1)./sqrt(nnz(issmallturn & hasonehs));

hs = [btd.hs];
acc = logical([hs.accepted]);
rej = ~acc;
isfirsths = [hs.hsnum] == 1;

% firstacc = acc(isfirsths);
% firstrej = rej(isfirsths);

num_hspts = 50;
dm = btd.behaviorTriggeredDataMatrix('hs', 'start', gqname, linspace(headSweepAcceptRange(1), headSweepAcceptRange(2), num_hspts));
dm2 = dm; dm2(~isfinite(dm2)) = 0;
hsmc = mean(dm2, 2); %rate per second
% firsthsmc = hsmc(isfirsths);

btdstruct.meanChangeAfterHeadsweep = mean(hsmc,1);
btdstruct.meanChangeAfterHeadsweep_eb = std(hsmc, 1)./sqrt(length(hsmc));
btdstruct.meanChangeAfterFirstHeadsweep = sum(hsmc(isfirsths,:), 1)./nnz(isfirsths);
btdstruct.meanChangeAfterFirstHeadsweep_eb = std(hsmc(isfirsths,:), 1)./sqrt(nnz(isfirsths));

btdstruct.meanChangeAfterAcceptedHeadsweep = sum(hsmc(acc,:), 1)./nnz(acc);
btdstruct.meanChangeAfterAcceptedHeadsweep_eb = std(hsmc(acc,:), 1)./sqrt(nnz(acc));

btdstruct.meanChangeAfterRejectedHeadsweep = sum(hsmc(rej,:), 1)./nnz(rej);
btdstruct.meanChangeAfterRejectedHeadsweep_eb = std(hsmc(rej,:), 1)./sqrt(nnz(rej));

btdstruct.meanChangeAfterFirstAcceptedHeadsweep = sum(hsmc(acc & isfirsths,:), 1)./nnz(acc & isfirsths);
btdstruct.meanChangeAfterFirstAcceptedHeadsweep_eb = std(hsmc(acc & isfirsths,:), 1)./sqrt(nnz(acc & isfirsths));

btdstruct.meanChangeAfterFirstRejectedHeadsweep = sum(hsmc(rej & isfirsths,:), 1)./nnz(rej & isfirsths);
btdstruct.meanChangeAfterFirstRejectedHeadsweep_eb = std(hsmc(rej & isfirsths,:), 1)./sqrt(nnz(rej & isfirsths));



firsths.name = 'hsnum';
firsths.validop = @(x) x == 1;
otherhs = firsths;
otherhs.validop = @(x) x > 1;

btdstruct.weightedFirstHeadsweep = btd.behaviorTriggeredWeightedAverage('hs', 'start', 'accepted', gqname, taxis, 'conditions', firsths);

btdstruct.btd = btd;
aa = [btd.all]; hsdt = median(diff([aa.eti]));
hstaxis = -10:hsdt:10;
hsbw = 0.5;
cumnsamp = 3;

if (~exist('IShannon_HShannon_initialization'))
    return;
end

[it.firsths.mi, it.firsths.tcent, it.firsths.bit_rate, it.firsths.tstart] = miFieldSignalVsTime (btd, 'hs', 'accepted', 'start', gqname, hstaxis, hsbw, 'conditions', firsths);
[~, ~, it.firsths.bit_rate_fast, it.firsths.tstart_cum] = miFieldSignalVsTime (btd, 'hs', 'accepted', 'start', gqname, hstaxis, [hstaxis(1:end-cumnsamp); hstaxis((1+cumnsamp):end)], 'conditions', firsths);
it.firsths.tend_cum = hstaxis((1+cumnsamp):end);
[it.otherhs.mi, it.otherhs.tcent, it.otherhs.bit_rate, it.otherhs.tstart] = miFieldSignalVsTime (btd, 'hs', 'accepted', 'start', gqname, hstaxis, hsbw, 'conditions', otherhs);
[~, ~, it.otherhs.bit_rate_fast, it.otherhs.tstart_cum] = miFieldSignalVsTime (btd, 'hs', 'accepted', 'start', gqname, hstaxis, [hstaxis(1:end-cumnsamp); hstaxis((1+cumnsamp):end)], 'conditions', otherhs);
it.otherhs.tend_cum = hstaxis((1+cumnsamp):end);
[it.allhs.mi, it.allhs.tcent, it.allhs.bit_rate, it.allhs.tstart] = miFieldSignalVsTime (btd, 'hs', 'accepted', 'start', gqname, hstaxis, hsbw);
[~, ~, it.allhs.bit_rate_fast, it.allhs.tstart_cum] = miFieldSignalVsTime (btd, 'hs', 'accepted', 'start', gqname, hstaxis, [hstaxis(1:end-cumnsamp); hstaxis((1+cumnsamp):end)]);
it.allhs.tend_cum = hstaxis((1+cumnsamp):end);
it.firsths.bw = hsbw;
it.allhs.bw = hsbw;
it.otherhs.bw = hsbw;


tstaxis = -60:0.25:60;
tsbw = 5;
cumnsamp = 2;
[it.turnsize.mi, it.turnsize.tcent, it.turnsize.bit_rate, it.turnsize.tstart] = miFieldSignalVsTime (btd, 'turn', 'dtheta', 'start', gqname, tstaxis, tsbw, 'operation', @(x) abs(x));
[~, ~, it.turnsize.bit_rate_fast, it.turnsize.tstart_cum] = miFieldSignalVsTime (btd, 'turn', 'dtheta', 'start', gqname, tstaxis, [tstaxis(1:end-cumnsamp); tstaxis((1+cumnsamp):end)], 'operation', @(x) abs(x));
it.turnsize.tend_cum = tstaxis((1+cumnsamp):end);
it.turnsize.bw = tsbw;
fn = {'firsths', 'otherhs', 'allhs','turnsize'};
for j = 1:length(fn)
    dt = diff(it.(fn{j}).tstart_cum); dt = dt([1:end end]);
    it.(fn{j}).cumbit = cumsum((it.(fn{j}).bit_rate_fast - msNoOutliers(it.(fn{j}).bit_rate_fast,2.5)).*dt);
    
    it.(fn{j}).totalbits = it.(fn{j}).cumbit(end);
    [it.(fn{j}).peakrate, I] = max((it.(fn{j}).bit_rate - msNoOutliers(it.(fn{j}).bit_rate,2.5)));
    it.(fn{j}).peaktime_start = it.(fn{j}).tstart(I);
    it.(fn{j}).peaktime_cent = it.(fn{j}).tcent(I);
    it.(fn{j}).peaktime_end = 2*it.(fn{j}).tcent(I) - it.(fn{j}).tstart(I);
    it.(fn{j}).mi = it.(fn{j}).mi/log(2); %provide mutual information in bits
end
it.unit = 'bits';
btdstruct.ShannonInfo = it;
    
