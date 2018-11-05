function btdstruct = analyzeBTDDirectory_twosource (btdstruct, stim1name, stim2name, varargin)
%function makeBTDFiguresBasic (basedir)

taxis = -30:0.1:10;
kernelTime = 7; %s
kernelDt = 0.1;
turnSizeRange = [-15 -3];
headSweepAcceptRange = [0 1.25];

vallist = {'led1Val', 'led2Val'};
varargin = assignApplicable(varargin);

if (~isstruct(btdstruct))
    if (ischar(btdstruct))
        btdstruct = BehaviorTriggeredData.loadBTDDirectory(btdstruct);
    else
        error ('first argument must be btdstruct or basedir');
    end
else
    if (isfield(btdstruct, 'stimnames'))
        existsAndDefault('stim1name', btdstruct.stimnames{1});
        existsAndDefault('stim2name', btdstruct.stimnames{2});
    end
    
end
% 
% turnSizeRange = [-15 -3];
% headSweepAcceptRange = [0 1.25];
btdstruct.turnSizeRange = turnSizeRange;
btdstruct.headSweepAcceptRange = headSweepAcceptRange;

btdstruct.kernelTime = kernelTime;
btdstruct.kernelDt = kernelDt;
btdstruct.taxis = taxis;

btd = btdstruct.btd;

kname1 = [stim1name 'TurnLin'];
kname2 = [stim2name 'TurnLin'];

inrun.name = 'isrun';
inrun.validop = @(x) logical(x);

[btdstruct.convkernels, btd] = btd.createBTAKernel('turn', 'start', {stim1name, stim2name}, kernelTime, kernelDt, 'newFieldName', {kname1, kname2}, 'normalizeStd', false, 'abbott', true);
aa = [btd.all];
eti = [aa.eti];
l1 = btd.behaviorTriggeredDataMatrix ('all', '', kname1, 0);
s = std(l1(eti > kernelTime));
btdstruct.convkernels(1,:) = btdstruct.convkernels(1,:) / s;
btd = btd.addConvolvedFields(stim1name, kname1, btdstruct.convkernels(1,:), kernelDt, 'scaleToSqr', true);

l2 = btd.behaviorTriggeredDataMatrix ('all', '', kname2, 0);
s = std(l2(eti > kernelTime));
btdstruct.convkernels(2,:) = btdstruct.convkernels(2,:) / s;
btd = btd.addConvolvedFields(stim2name, kname2, btdstruct.convkernels(2,:), kernelDt, 'scaleToSqr', true);

xyall =  btd.behaviorTriggeredDataMatrix('all', 'start', {kname1, kname2}, 0);
mx = mean(xyall(eti>kernelTime,1)); my = mean(xyall(eti>kernelTime,1));
sx = std(xyall(eti>kernelTime,2)); sy = std(xyall(eti>kernelTime,2));

lin_xy_all = btd.behaviorTriggeredDataMatrix('all', 'start', {kname1, kname2}, 0,'conditions', inrun);
lin_xy_all(:,1) = (lin_xy_all(:,1)-mx)/sx;
lin_xy_all(:,2) = (lin_xy_all(:,2)-my)/sy;


btdstruct.lin_xy = btd.behaviorTriggeredDataMatrix('turn', 'start', {kname1, kname2}, 0);
btdstruct.lin_xy_all = lin_xy_all;

theta = atan2(mean(btdstruct.lin_xy(:,2)), mean(btdstruct.lin_xy(:,1)));

btdstruct.theta = theta;

btdstruct.linUOp =  @(yd) cos(theta)*(yd{1}-mx)/sx + sin(theta)*(yd{2}-my)/sy;
btdstruct.linVOp = @(yd) cos(theta)*(yd{2}-my)/sy - sin(theta)*(yd{1}-mx)/sx;
% btd = addOperationFields (btd, gltname, newname, operation)
btd = btd.addOperationFields({kname1, kname2}, 'linU',btdstruct.linUOp);
btd = btd.addOperationFields({kname1, kname2}, 'linV', btdstruct.linVOp );

m = sign(sum(btdstruct.convkernels, 2));

btd = btd.addOperationFields({stim1name, stim2name}, 'ufieldsum', @(yd) m(1) * cos(theta)*yd{1} + m(2)*sin(theta)*yd{2});
btd = btd.addOperationFields({stim1name, stim2name}, 'vfieldsum', @(yd) m(1) *cos(theta)*yd{2} - m(2)*sin(theta)*yd{1});


btdstruct.lin_uv = btd.behaviorTriggeredDataMatrix('turn', 'start', {'linU', 'linV'}, 0);
%btdstruct.lin_uv_all = btd.behaviorTriggeredDataMatrix('all', 'start', {'linU', 'linV'}, 0, 'conditions', inrun);
% scale = repmat(std(lin_xy_all), [length(btdstruct.lin_xy) 1]);
% lin_xy = btdstruct.lin_xy./scale;

lin_xy(:,1) = (btdstruct.lin_xy(:,1)-mx)/sx;
lin_xy(:,2) = (btdstruct.lin_xy(:,2)-my)/sy;
lin_uv = btdstruct.lin_uv;


u = mean(lin_xy);
c = cov(lin_xy);
btdstruct.kl_xy = 0.5*(trace(c)-log(det(c))+u*u' - 2);

c = cov(lin_xy(:,1));
btdstruct.kl_x = 0.5*(trace(c)-log(det(c))+u(1).^2 - 1);

c = cov(lin_xy(:,2));
btdstruct.kl_y = 0.5*(trace(c)-log(det(c))+u(2).^2 - 1);


u = mean(lin_uv);
c = cov(lin_uv);
btdstruct.kl_uv = 0.5*(trace(c)-log(det(c))+u*u' - 2);

c = cov(lin_uv(:,1));
btdstruct.kl_u = 0.5*(trace(c)-log(det(c))+u(1).^2 - 1);

c = cov(lin_uv(:,2));
btdstruct.kl_v = 0.5*(trace(c)-log(det(c))+u(2).^2 - 1);
%find theta by minimizing u*v
minfun = @(t) (mean((cos(t)*lin_xy(:,1) + sin(t)*lin_xy(:,2)).*(-sin(t)*lin_xy(:,1) + cos(t)*lin_xy(:,2))))^2;
theta = fminunc(minfun, btdstruct.theta);

btdstruct.thetaAlt = theta;

btdstruct.linUOpAlt =  @(yd) cos(theta)*(yd{1}-mx)/sx + sin(theta)*(yd{2}-my)/sy;
btdstruct.linVOpAlt = @(yd) cos(theta)*(yd{2}-my)/sy - sin(theta)*(yd{1}-mx)/sx;
% btd = addOperationFields (btd, gltname, newname, operation)
btd = btd.addOperationFields({kname1, kname2}, 'linUAlt',btdstruct.linUOpAlt);
btd = btd.addOperationFields({kname1, kname2}, 'linVAlt', btdstruct.linVOpAlt );

% btd = btd.addOperationFields({stim1name, stim2name}, 'ufieldsum', @(yd) cos(theta)*yd{1} + sin(theta)*yd{2});
% btd = btd.addOperationFields({stim1name, stim2name}, 'vfieldsum', @(yd) cos(theta)*yd{2} - sin(theta)*yd{1});


btdstruct.lin_uvAlt = btd.behaviorTriggeredDataMatrix('turn', 'start', {'linUAlt', 'linVAlt'}, 0);
%btdstruct.lin_uv_all = btd.behaviorTriggeredDataMatrix('all', 'start', {'linU', 'linV'}, 0, 'conditions', inrun);
% scale = repmat(std(lin_xy_all), [length(btdstruct.lin_xy) 1]);
% lin_xy = btdstruct.lin_xy./scale;

lin_uv = btdstruct.lin_uvAlt;


u = mean(lin_uv);
c = cov(lin_uv);
btdstruct.kl_uvAlt = 0.5*(trace(c)-log(det(c))+u*u' - 2);

c = cov(lin_uv(:,1));
btdstruct.kl_uAlt = 0.5*(trace(c)-log(det(c))+u(1).^2 - 1);

c = cov(lin_uv(:,2));
btdstruct.kl_vAlt = 0.5*(trace(c)-log(det(c))+u(2).^2 - 1);



fieldlist = {kname1, kname2, 'linU', 'linV'};
ratenames = {'rateX', 'rateY', 'rateU', 'rateV'};
%cumnames = {'cumX', 'cumY', 'cumU', 'cumV'};
turndata = btd.behaviorTriggeredDataMatrix('turn', 'start', fieldlist, 0);
alldata =  btd.behaviorTriggeredDataMatrix('all', 'start', fieldlist, 0, 'conditions', inrun);
nr = size(turndata, 1);
na = size(alldata,1);
mu = mean(turndata);
mu_a = mean(alldata);
s = std(turndata);
s_a = std(alldata);
dt = median(diff(btd(1).all.eti)); %assuming all movies were captured at same frame rate!

for j = 1:size(turndata, 2)
    ratefun = @(ydata) nr/(na*dt)*normpdf(ydata,mu(j),s(j))./normpdf(ydata, mu_a(j), s_a(j));
    btdstruct.ratefun{j} = ratefun;
    btd = btd.addOperationFields(fieldlist{j}, ratenames{j}, ratefun);
%     
%     for k = 1:length(btd)
%         glt = btd(k).glt(btd(k).findField(ratenames{j}));
%         glt.fieldname = cumnames{j};
%         glt.yData = cumsum([0 diff(glt.xData)].*glt.yData);
%         ind = btd(k).findField(cumnames{j});
%         if (ind <= 0)
%             btd(k).glt = [btd(k).glt glt];
%         else
%             btd(k).glt(ind) = glt;
%         end
%     end
end


inrun.name = 'isrun';
inrun.validop = @(x) logical(x);


allxy = btd.behaviorTriggeredDataMatrix('all', '', {kname1, kname2}, 0, 'conditions', inrun);
turnxy = btd.behaviorTriggeredDataMatrix('turn', 'start',  {kname1, kname2}, 0);

dt = median(diff(btd(1).all.eti));

ratefun1 = @(x,xdata) x(1) * normpdf(xdata(:,1), x(2), x(3))./normpdf(xdata(:,1), 0, 1) + x(4)*normpdf(xdata(:,2), x(5), x(6))./normpdf(xdata(:,2), 0, 1);
initguess1 = [0.5*length(turnxy)/length(allxy) mean(turnxy(:,1)) std(turnxy(:,1)) 0.5*length(turnxy)/length(allxy) mean(turnxy(:,2)) std(turnxy(:,2))];

ratefun2 = @(x, xdata) x(1) * normpdf(cos(x(4))*xdata(:,1) + sin(x(4))*xdata(:,2), x(2), x(3)) ./ normpdf(cos(x(4))*xdata(:,1) + sin(x(4))*xdata(:,2), 0, 1); %fit rate, not theta
initguess2 = [0.5*length(turnxy)/length(allxy) mean(lin_uv(:,1)) std(lin_uv(:,1)) btdstruct.theta];

ratefun3 = @(x,xdata) x(1) * ( normpdf(xdata(:,1), x(2), x(3))./normpdf(xdata(:,1), 0, 1) ) .* ( normpdf(xdata(:,2), x(4), x(5))./normpdf(xdata(:,2), 0, 1) );
initguess3 = [0.5*length(turnxy)/length(allxy) mean(turnxy(:,1)) std(turnxy(:,1)) mean(turnxy(:,2)) std(turnxy(:,2))];



%if we assume rdt << 1
fitfun1 = @(x) -sum(log(max(ratefun1(x,turnxy),1E-100))) + sum(max(0,ratefun1(x,allxy)));
fitfun2 = @(x) -sum(log(max(ratefun2(x,turnxy),1E-100))) + sum(max(0,ratefun2(x,allxy)));
fitfun3 = @(x) -sum(log(max(ratefun3(x,turnxy),1E-100))) + sum(max(0,ratefun3(x,allxy)));

op = optimset('fmincon');
op.Algorithm = 'active-set';
%op.LargeScale = 'off';
problem.objective = fitfun2;
problem.solver = 'fmincon';
problem.x0 = initguess2;
problem.lb = [0 -Inf 0 0];
problem.ub = [Inf Inf Inf pi];
problem.options = op;
fitParams2 = fmincon(problem);
btdstruct.ratefun2 = ratefun2;
btdstruct.fitParams2 = fitParams2 .* [1/dt 1 1 1];
btdstruct.llAll2 = -fitfun2(fitParams2);

op = optimset('fmincon');
op.Algorithm = 'active-set';
%op.LargeScale = 'off';
problem.objective = fitfun1;
problem.solver = 'fmincon';
problem.x0 = initguess1;
problem.lb = [0 -Inf 0 0 -Inf 0];
problem.ub = [Inf Inf Inf Inf Inf Inf];
problem.options = op;
fitParams1 = fmincon(problem);
btdstruct.ratefun1 = ratefun1;
btdstruct.fitParams1 = fitParams1.*[1/dt 1 1 1/dt 1 1];
btdstruct.llAll1 = -fitfun1(fitParams1);

op = optimset('fmincon');
op.Algorithm = 'active-set';
%op.LargeScale = 'off';
problem.objective = fitfun3;
problem.solver = 'fmincon';
problem.x0 = initguess3;
problem.lb = [0 -Inf 0 -Inf 0];
problem.ub = [Inf Inf Inf Inf Inf];
problem.options = op;
fitParams3 = fmincon(problem);
btdstruct.ratefun3 = ratefun3;
btdstruct.fitParams3 = fitParams3.*[1/dt 1 1 1 1];
btdstruct.llAll3 = -fitfun3(fitParams3);
    
    

btdstruct.btd = btd;
btdstruct.fieldlist = fieldlist;
btdstruct.ratenames = ratenames;
% btdstruct.cumnames = cumnames;
% 
% dm = btd.behaviorTriggeredDataMatrix('turn', 'start', cumnames, 0);
% tt = [btd.turn]; tn = [tt.tnum];    
% valid = (diff(tn) == 0);
% dm = diff(dm);
% % 
% for j = 1:length(cumnames)
%     btdstruct.lambda(:,j) = dm(valid,j);
% end


%link = btdstruct.btd.behaviorTriggeredDataMatrix('turn', 'start', {kname1, kname2}, 0);

firsths.name = 'hsnum';
firsths.validop = @(x) x == 1;


quadrant(1).condition(1).name = kname1;
quadrant(1).condition(1).validop = @(x) x > 0;
quadrant(1).condition(2).name = kname2;
quadrant(1).condition(2).validop = @(x) x > 0;

%changed quadrants 10/28 to align with normal cartesian labels
quadrant(2) = quadrant(1);
quadrant(2).condition(1).validop = @(x) x < 0;

quadrant(3) = quadrant(2);
quadrant(3).condition(2).validop = @(x) x < 0;

quadrant(4) = quadrant(1);
quadrant(4).condition(2).validop = @(x) x < 0;


plane(1).condition.name = kname1;
plane(1).condition.validop = @(x) x>0;
plane(2).condition.name = kname1;
plane(2).condition.validop = @(x) x<0;
plane(3:4) = plane(1:2);
plane(3).condition.name = kname2;
plane(4).condition.name = kname2;

btdstruct.taxis = taxis;


btdstruct.stimnames = {stim1name stim2name};
btdstruct.linnames = {kname1, kname2};

t = [btd.turn]; 
btdstruct.meanSquaredTurn = mean([t.dtheta].^2);
ms = btdstruct.meanSquaredTurn;
bigturn.name = 'dtheta';
bigturn.validop = @(x) x.^2 > ms;
smallturn.name = 'dtheta';
smallturn.validop = @(x) x.^2 < ms;
onehs.name = 'numhs';
onehs.validop = @(x) x == 1;
stimlist = {stim1name,stim2name, 'ufieldsum', 'vfieldsum'};%,'linU', 'linV',kname1, kname2};
wholelist = [stimlist vallist];

num_tspts = 200;

[dm,stimnum] = btd.behaviorTriggeredDataMatrix('turn', 'start', stimlist, linspace(turnSizeRange(1), turnSizeRange(2), num_tspts));
turnmc = zeros(size(dm,1), length(stimlist));

for j = 1:length(stimlist)
    dm2 = reshape(dm(stimnum == j), size(dm,1),[]); dm2(~isfinite(dm2)) = 0;
    turnmc(:,j) = mean(dm2, 2); %rate per second
end
t = [btd.turn]; 
isbigturn = ([t.dtheta].^2 > btdstruct.meanSquaredTurn);
issmallturn = ~isbigturn;
hasonehs = [t.numhs] == 1;

btdstruct.all.meanChangeBeforeTurn = mean(turnmc, 1);
btdstruct.all.meanChangeBeforeTurn_eb = std(turnmc, 1)./sqrt(length(turnmc));

btdstruct.all.meanChangeBeforeTurn_onehs = sum(turnmc(hasonehs,:), 1)./nnz(hasonehs);
btdstruct.all.meanChangeBeforeTurn_onehs_eb = std(turnmc(hasonehs,:), 1)./sqrt(nnz(hasonehs));

btdstruct.all.meanChangeBeforeBigTurn = sum(turnmc(isbigturn,:), 1)./nnz(isbigturn);
btdstruct.all.meanChangeBeforeBigTurn_eb = std(turnmc(isbigturn,:), 1)./sqrt(nnz(isbigturn));

btdstruct.all.meanChangeBeforeSmallTurn = sum(turnmc(issmallturn,:), 1)./nnz(issmallturn);
btdstruct.all.meanChangeBeforeSmallTurn_eb = std(turnmc(issmallturn,:), 1)./sqrt(nnz(issmallturn));

btdstruct.all.meanChangeBeforeBigTurn_onehs = sum(turnmc(isbigturn & hasonehs,:), 1)./nnz(isbigturn & hasonehs);
btdstruct.all.meanChangeBeforeBigTurn_onehs_eb = std(turnmc(isbigturn & hasonehs,:), 1)./sqrt(nnz(isbigturn & hasonehs));

btdstruct.all.meanChangeBeforeSmallTurn_onehs = sum(turnmc(issmallturn & hasonehs,:), 1)./nnz(issmallturn & hasonehs);
btdstruct.all.meanChangeBeforeSmallTurn_onehs_eb = std(turnmc(issmallturn & hasonehs,:), 1)./sqrt(nnz(issmallturn & hasonehs));

hs = [btd.hs];
acc = logical([hs.accepted]);
rej = ~acc;
isfirsths = [hs.hsnum] == 1;

firstacc = acc(isfirsths);
firstrej = rej(isfirsths);

num_hspts = 40;
[dm,stimnum] = btd.behaviorTriggeredDataMatrix('hs', 'start', stimlist, linspace(headSweepAcceptRange(1), headSweepAcceptRange(2), num_hspts));
hsmc = zeros(size(dm,1), length(stimlist));
firsthsmc = zeros(nnz(isfirsths), length(stimlist));
for j = 1:length(stimlist)
    dm2 = reshape(dm(stimnum == j), size(dm,1),[]); dm2(~isfinite(dm2)) = 0;
    hsmc(:,j) = mean(dm2, 2); %rate per second
    firsthsmc(:,j) = hsmc(isfirsths, j);
end

btdstruct.all.meanChangeAfterHeadsweep = mean(hsmc,1);
btdstruct.all.meanChangeAfterHeadsweep_eb = std(hsmc, 1)./sqrt(length(hsmc));
btdstruct.all.meanChangeAfterFirstHeadsweep = sum(hsmc(isfirsths,:), 1)./nnz(isfirsths);
btdstruct.all.meanChangeAfterFirstHeadsweep_eb = std(hsmc(isfirsths,:), 1)./sqrt(nnz(isfirsths));


btdstruct.all.meanChangeAfterAcceptedHeadsweep = sum(hsmc(acc,:), 1)./nnz(acc);
btdstruct.all.meanChangeAfterAcceptedHeadsweep_eb = std(hsmc(acc,:), 1)./sqrt(nnz(acc));

btdstruct.all.meanChangeAfterRejectedHeadsweep = sum(hsmc(rej,:), 1)./nnz(rej);
btdstruct.all.meanChangeAfterRejectedHeadsweep_eb = std(hsmc(rej,:), 1)./sqrt(nnz(rej));

btdstruct.all.meanChangeAfterFirstAcceptedHeadsweep = sum(hsmc(acc & isfirsths,:), 1)./nnz(acc & isfirsths);
btdstruct.all.meanChangeAfterFirstAcceptedHeadsweep_eb = std(hsmc(acc & isfirsths,:), 1)./sqrt(nnz(acc & isfirsths));

btdstruct.all.meanChangeAfterFirstRejectedHeadsweep = sum(hsmc(rej & isfirsths,:), 1)./nnz(rej & isfirsths);
btdstruct.all.meanChangeAfterFirstRejectedHeadsweep_eb = std(hsmc(rej & isfirsths,:), 1)./sqrt(nnz(rej & isfirsths));

fields = {'turn', 'acchs', 'rejhs'};% 'acchs', 'rejhs','acchs', 'rejhs'};
pos =   {'start', 'start', 'start'};% 'atmax', 'atmax',  'end',   'end'};

for j = 1:length(fields)
    [bta, ~] = btd.behaviorTriggeredAverage(fields{j}, pos{j},wholelist, taxis);       
    btdstruct.all.([fields{j} '_' pos{j}]) = reshape(bta, [], length(wholelist));
end
[bta,~,eb] = btd.behaviorTriggeredAverage('turn', 'start',wholelist, taxis, 'conditions', [bigturn onehs]);
btdstruct.all.bigturn_start_onehs = reshape(bta,[],length(wholelist));
btdstruct.all.bigturn_start_onehs_eb = reshape(eb,[],length(wholelist));
[bta,~,eb] = btd.behaviorTriggeredAverage('turn', 'start',wholelist, taxis, 'conditions', [smallturn onehs]);
btdstruct.all.smallturn_start_onehs = reshape(bta,[],length(wholelist));
btdstruct.all.smallturn_start_onehs_eb = reshape(eb,[],length(wholelist));
[bta,~,eb] = btd.behaviorTriggeredAverage('turn', 'start',wholelist, taxis, 'conditions', bigturn );
btdstruct.all.bigturn_start = reshape(bta,[],length(wholelist));
btdstruct.all.bigturn_start_eb = reshape(eb,[],length(wholelist));
[bta,~,eb] = btd.behaviorTriggeredAverage('turn', 'start',wholelist, taxis, 'conditions', smallturn );
btdstruct.all.smallturn_start = reshape(bta,[],length(wholelist));
btdstruct.all.smallturn_start_eb = reshape(eb,[],length(wholelist));

fields = {'turn'};% 'acchs', 'rejhs','acchs', 'rejhs'};
pos =   {'start'};% 'atmax', 'atmax',  'end',   'end'};

for k = 1:length(btd)
    
    for j = 1:length(fields)
        btdstruct.exp(k).([fields{j} '_' pos{j}]) = reshape(btd(k).behaviorTriggeredAverage(fields{j}, pos{j}, wholelist, taxis),[],length(wholelist));
    end
end
% 
% btdstruct.all.smallturn_start_onehs = reshape(btd.behaviorTriggeredAverage('turn', 'start',wholelist, taxis, 'conditions', [smallturn onehs]),[],length(wholelist));
% btdstruct.all.bigturn_start = reshape(btd.behaviorTriggeredAverage('turn', 'start',wholelist, taxis, 'conditions', bigturn),[],length(wholelist));
% btdstruct.all.smallturn_start = reshape(btd.behaviorTriggeredAverage('turn', 'start',wholelist, taxis, 'conditions', smallturn),[],length(wholelist));


fields = {'acchs', 'rejhs'};%, 'acchs', 'rejhs','acchs', 'rejhs'};
pos =   {'start', 'start'};%, 'atmax', 'atmax',  'end',   'end'};

% function [bta, dm] = behaviorTriggeredAverage (btd, tp, position, gltname, displacementAxis, varargin)
for j = 1:length(fields)
    [bta, ~] = btd.behaviorTriggeredAverage(fields{j}, pos{j},wholelist, taxis, 'conditions', firsths);
    %btdstruct.([fields{j} '_' pos{j} '_dm']) = dm; %this takes up a lot of space!
    btdstruct.all.(['first_' fields{j} '_' pos{j}]) = reshape(bta, [], length(wholelist));
end


btdstruct.all.weightedTurnStart_onehs = reshape(btd.behaviorTriggeredWeightedAverage('turn', 'start', 'dtheta', wholelist, taxis, 'weightedOp', @abs, 'conditions', onehs), [], length(wholelist));
btdstruct.all.weightedTurnStart_all = reshape(btd.behaviorTriggeredWeightedAverage('turn', 'start', 'dtheta', wholelist, taxis, 'weightedOp', @abs), [], length(wholelist));

firsths.name = 'hsnum';
firsths.validop = @(x) x == 1;
btdstruct.all.weightedFirstHeadsweep = reshape(btd.behaviorTriggeredWeightedAverage('hs', 'start', 'accepted', wholelist, taxis, 'conditions', firsths), [], length(wholelist));


kv = btd.behaviorTriggeredDataMatrix('turn', 'start', {quadrant(1).condition.name}, 0);

% for j = 1:length(quadrant(1).condition)
%     kv(:,j) = btd.behaviorTriggeredDataMatrix('turn', 'start', quadrant(1).condition(j).name, 0);
% end
qq = zeros(size(kv,1), 1);
for k = 1:length(quadrant)
    valid = true(size(qq));
    for j = 1:length(quadrant(k).condition)
        valid = valid & quadrant(k).condition(j).validop(kv(:,j));
    end
    if (any(nnz(qq(valid))))
        warning('same point satisfies multiple quadrants'); %#ok<WNTAG>
    end
    qq(valid) = k;
end
    
qq = qq';

for k = 1:length(quadrant)
    

    fields = {'turn', 'acchs', 'rejhs'};% 'acchs', 'rejhs','acchs', 'rejhs'};
    pos =   {'start', 'start', 'start'};% 'atmax', 'atmax',  'end',   'end'};


    for j = 1:length(fields)
        [bta, ~, eb] = btd.behaviorTriggeredAverage(fields{j}, pos{j}, wholelist, taxis, 'conditions', quadrant(k).condition);       
        quadrant(k).([fields{j} '_' pos{j}]) = reshape(bta, [], length(wholelist));
        quadrant(k).([fields{j} '_' pos{j} '_eb']) = reshape(eb, [], length(wholelist));
    end
    
    fields = {'acchs', 'rejhs'};%, 'acchs', 'rejhs','acchs', 'rejhs'};
    pos =   {'start', 'start'};%, 'atmax', 'atmax',  'end',   'end'};

    % function [bta, dm] = behaviorTriggeredAverage (btd, tp, position, gltname, displacementAxis, varargin)
    for j = 1:length(fields)
        [bta, ~, eb] = btd.behaviorTriggeredAverage(fields{j}, pos{j}, wholelist, taxis, 'conditions', [quadrant(k).condition firsths]);
        %btdstruct.([fields{j} '_' pos{j} '_dm']) = dm; %this takes up a lot of space!
        quadrant(k).(['first_' fields{j} '_' pos{j}]) = reshape(bta, [], length(wholelist));
        quadrant(k).(['first_' fields{j} '_' pos{j} '_eb']) = reshape(eb, [], length(wholelist));
    end

    [bta,~,eb] = btd.behaviorTriggeredAverage('turn', 'start',wholelist, taxis, 'conditions', [quadrant(k).condition bigturn onehs]);
    quadrant(k).bigturn_start_onehs = reshape(bta,[],length(wholelist));
    quadrant(k).bigturn_start_onehs_eb = reshape(eb,[],length(wholelist));
    [bta,~,eb] = btd.behaviorTriggeredAverage('turn', 'start',wholelist, taxis, 'conditions', [quadrant(k).condition smallturn onehs]);
   quadrant(k).smallturn_start_onehs = reshape(bta,[],length(wholelist));
    quadrant(k).smallturn_start_onehs_eb = reshape(eb,[],length(wholelist));
    [bta,~,eb] = btd.behaviorTriggeredAverage('turn', 'start',wholelist, taxis, 'conditions', [quadrant(k).condition bigturn] );
    quadrant(k).bigturn_start = reshape(bta,[],length(wholelist));
    quadrant(k).bigturn_start_eb = reshape(eb,[],length(wholelist));
    [bta,~,eb] = btd.behaviorTriggeredAverage('turn', 'start',wholelist, taxis, 'conditions', [quadrant(k).condition smallturn] );
    quadrant(k).smallturn_start = reshape(bta,[],length(wholelist));
    quadrant(k).smallturn_start_eb = reshape(eb,[],length(wholelist));
%     quadrant(k).bigturn_start_onehs = reshape(btd.behaviorTriggeredAverage('turn', 'start', wholelist, taxis, 'conditions', [quadrant(k).condition bigturn onehs]),[],length(wholelist));
%     quadrant(k).smallturn_start_onehs = reshape(btd.behaviorTriggeredAverage('turn', 'start', wholelist, taxis, 'conditions', [quadrant(k).condition smallturn onehs]),[],length(wholelist));
%     quadrant(k).weightedTurnStart_onehs =  reshape(btd.behaviorTriggeredWeightedAverage('turn', 'start', 'dtheta', wholelist, taxis, 'weightedOp', @abs, 'conditions', [quadrant(k).condition onehs]),[],length(wholelist));
    quadrant(k).weightedFirstHeadsweep =  reshape(btd.behaviorTriggeredWeightedAverage('hs', 'start', 'accepted', wholelist, taxis, 'conditions',  [quadrant(k).condition firsths]),[],length(wholelist));
    
    quadrant(k).meanChangeBeforeTurn = sum(turnmc(qq == k,:), 1)./nnz(qq == k);
    quadrant(k).meanChangeBeforeTurn_eb = std(turnmc(qq == k,:), 1)./sqrt(nnz(qq == k));
    quadrant(k).meanChangeBeforeTurn_onehs = sum(turnmc(hasonehs & qq == k,:), 1)./nnz(hasonehs & qq == k);
    quadrant(k).meanChangeBeforeTurn_onehs_eb = std(turnmc(hasonehs & qq == k,:), 1)./sqrt(nnz(hasonehs & qq == k));

    
    quadrant(k).meanChangeBeforeBigTurn = sum(turnmc(isbigturn & qq == k,:), 1)./nnz(isbigturn & qq == k);
    quadrant(k).meanChangeBeforeBigTurn_eb = std(turnmc(isbigturn & qq == k,:), 1)./sqrt(nnz(isbigturn & qq == k));

    quadrant(k).meanChangeBeforeSmallTurn = sum(turnmc(issmallturn & qq == k,:), 1)./nnz(issmallturn & qq == k);
    quadrant(k).meanChangeBeforeSmallTurn_eb = std(turnmc(issmallturn & qq == k,:), 1)./sqrt(nnz(issmallturn & qq == k));

    quadrant(k).meanChangeBeforeBigTurn_onehs = sum(turnmc(isbigturn & hasonehs & qq == k,:), 1)./nnz(isbigturn & hasonehs & qq == k);
    quadrant(k).meanChangeBeforeBigTurn_onehs_eb = std(turnmc(isbigturn & hasonehs & qq == k,:), 1)./sqrt(nnz(isbigturn & hasonehs & qq == k));

    quadrant(k).meanChangeBeforeSmallTurn_onehs = sum(turnmc(issmallturn & hasonehs & qq == k,:), 1)./nnz(issmallturn & hasonehs & qq == k);
    quadrant(k).meanChangeBeforeSmallTurn_onehs_eb = std(turnmc(issmallturn & hasonehs & qq == k,:), 1)./sqrt(nnz(issmallturn & hasonehs & qq == k));

    quadrant(k).meanChangeAfterFirstHeadsweep = sum(firsthsmc(qq == k,:), 1)./nnz(qq == k);
    quadrant(k).meanChangeAfterFirstHeadsweep_eb = std(firsthsmc(qq == k,:), 1)./sqrt(nnz(qq == k));
     
    quadrant(k).meanChangeAfterFirstAcceptedHeadsweep = sum(firsthsmc(firstacc & qq == k,:), 1)./nnz(firstacc & qq == k);
    quadrant(k).meanChangeAfterFirstAcceptedHeadsweep_eb = std(firsthsmc(firstacc & qq == k,:), 1)./sqrt(nnz(firstacc & qq == k));

    quadrant(k).meanChangeAfterFirstRejectedHeadsweep = sum(firsthsmc(firstrej & qq == k,:), 1)./nnz(firstrej & qq == k);
    quadrant(k).meanChangeAfterFirstRejectedHeadsweep_eb = std(firsthsmc(firstrej & qq == k,:), 1)./sqrt(nnz(firstrej & qq == k));

end

btdstruct.quadrant = quadrant;

return;
for k = 1:length(plane)
    

    fields = {'turn', 'acchs', 'rejhs'};% 'acchs', 'rejhs','acchs', 'rejhs'};
    pos =   {'start', 'start', 'start'};% 'atmax', 'atmax',  'end',   'end'};


    for j = 1:length(fields)
        [bta, ~,eb] = btd.behaviorTriggeredAverage(fields{j}, pos{j}, wholelist, taxis, 'conditions', plane(k).condition);       
        plane(k).([fields{j} '_' pos{j}]) = reshape(bta, [], length(wholelist));
        plane(k).([fields{j} '_' pos{j} '_eb']) = reshape(eb, [], length(wholelist));
    end
    
    fields = {'acchs', 'rejhs'};%, 'acchs', 'rejhs','acchs', 'rejhs'};
    pos =   {'start', 'start'};%, 'atmax', 'atmax',  'end',   'end'};

    % function [bta, dm] = behaviorTriggeredAverage (btd, tp, position, gltname, displacementAxis, varargin)
    for j = 1:length(fields)
        [bta, ~,eb] = btd.behaviorTriggeredAverage(fields{j}, pos{j}, wholelist, taxis, 'conditions', [plane(k).condition firsths]);
        %btdstruct.([fields{j} '_' pos{j} '_dm']) = dm; %this takes up a lot of space!
        plane(k).(['first_' fields{j} '_' pos{j}]) = reshape(bta, [], length(wholelist));
        plane(k).(['first_' fields{j} '_' pos{j} '_eb']) = reshape(eb, [], length(wholelist));
        
    end

    
    plane(k).bigturn_start_onehs = reshape(btd.behaviorTriggeredAverage('turn', 'start', wholelist, taxis, 'conditions', [plane(k).condition bigturn onehs]),[],length(wholelist));
    plane(k).smallturn_start_onehs = reshape(btd.behaviorTriggeredAverage('turn', 'start', wholelist, taxis, 'conditions', [plane(k).condition smallturn onehs]),[],length(wholelist));
    plane(k).weightedTurnStart_onehs = reshape(btd.behaviorTriggeredWeightedAverage('turn', 'start', 'dtheta', wholelist, taxis, 'weightedOp', @abs, 'conditions', [plane(k).condition onehs]),[],length(wholelist));
    plane(k).weightedFirstHeadsweep = reshape(btd.behaviorTriggeredWeightedAverage('hs', 'start', 'accepted', wholelist, taxis, 'conditions',  [plane(k).condition firsths]),[],length(wholelist));

end

btdstruct.plane = plane;

%

% 

% 
% 
% theta = atan2(mean(btdstruct.lin_xy(:,2)), mean(btdstruct.lin_xy(:,1)));
% 
% btdstruct.theta_old = theta;
% 
% btdstruct.linUOpOld =  @(yd) cos(theta)*(yd{1}-mx)/sx + sin(theta)*(yd{2}-my)/sy;
% btdstruct.linVOpOld = @(yd) cos(theta)*(yd{2}-my)/sy - sin(theta)*(yd{1}-mx)/sx;
% % btd = addOperationFields (btd, gltname, newname, operation)
% btd = btd.addOperationFields({kname1, kname2}, 'linUOld',btdstruct.linUOpOld);
% btd = btd.addOperationFields({kname1, kname2}, 'linVOld', btdstruct.linVOpOld);
% 
% % btd = btd.addOperationFields({stim1name, stim2name}, 'ufieldsum', @(yd) cos(theta)*yd{1} + sin(theta)*yd{2});
% % btd = btd.addOperationFields({stim1name, stim2name}, 'vfieldsum', @(yd) cos(theta)*yd{2} - sin(theta)*yd{1});
% 
% 
% btdstruct.lin_uv_old = btd.behaviorTriggeredDataMatrix('turn', 'start', {'linUOld', 'linVOld'}, 0);
% %btdstruct.lin_uv_all = btd.behaviorTriggeredDataMatrix('all', 'start', {'linU', 'linV'}, 0, 'conditions', inrun);
% % scale = repmat(std(lin_xy_all), [length(btdstruct.lin_xy) 1]);
% % lin_xy = btdstruct.lin_xy./scale;
% 
% lin_xy(:,1) = (btdstruct.lin_xy(:,1)-mx)/sx;
% lin_xy(:,2) = (btdstruct.lin_xy(:,2)-my)/sy;
% lin_uv = btdstruct.lin_uv_old;
% 
% 
% u = mean(lin_xy);
% c = cov(lin_xy);
% btdstruct.kl_xy = 0.5*(trace(c)-log(det(c))+u*u' - 2);
% 
% c = cov(lin_xy(:,1));
% btdstruct.kl_x = 0.5*(trace(c)-log(det(c))+u(1).^2 - 1);
% 
% c = cov(lin_xy(:,2));
% btdstruct.kl_y = 0.5*(trace(c)-log(det(c))+u(2).^2 - 1);
% 
% 
% u = mean(lin_uv);
% c = cov(lin_uv);
% btdstruct.kl_uv_old = 0.5*(trace(c)-log(det(c))+u*u' - 2);
% 
% c = cov(lin_uv(:,1));
% btdstruct.kl_u_old = 0.5*(trace(c)-log(det(c))+u(1).^2 - 1);
% 
% c = cov(lin_uv(:,2));
% btdstruct.kl_v_old = 0.5*(trace(c)-log(det(c))+u(2).^2 - 1);