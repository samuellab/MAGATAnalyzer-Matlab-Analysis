function btdstruct = prepVarianceSwitchingAnalysis(btdstruct, opstruct)
%btdstruct = prepVarianceSwitchingAnalysis(btdstruct, opstruct)
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
    btdstruct = opstruct;
    disp ('ramp types are square, triangle, sine, and constant');
    return;
end

if (nargin < 2 || ~isfield(btdstruct, 'btd'))
    warning ('first argument must be a btdstruct and second argument contains options');
    return;
end


for j = find([opstruct.stim.iswn])
    if (opstruct.stim(j).switch == 1)
        btdstruct.btd = btdstruct.btd.addVarianceGQs(opstruct.stim(j).gqname, opstruct.stim(j).ramptype, 'fixedPeriod', opstruct.stim(j).period);
    end
   % if (opstruct.redogq || any(btdstruct.btd.findField([opstruct.stim(j).gqname '_var_ton']) <= 0))
    %end
    
end

for j = find(~[opstruct.stim.iswn] & [opstruct.stim.switch])
    %if (opstruct.redogq || any(btdstruct.btd.findField([opstruct.stim(j).gqname '_ton']) <= 0))
        btdstruct.btd = btdstruct.btd.addTonToffGQs(opstruct.stim(j).gqname, opstruct.stim(j).ramptype, 'fixedPeriod', opstruct.stim(j).period);
    %end
end

for j = find(~[opstruct.stim.iswn] & [opstruct.stim.switch])
    disp (['stim ' - num2str(j) ' is not switching or noisy - so I am ignoring it']);
end

varswitch = [opstruct.stim.iswn] & [opstruct.stim.switch];
switchstim = 0;
switch (nnz(varswitch))
    case 0
        if (any([opstruct.stim.switch]))
            switchstim = find([opstruct.stim.switch], 1, 'first');
            gq_switch = opstruct.stim(switchstim).gqname;
        end
    case 1
%         gq_switch = opstruct.stim(varswitch).gqname;
        switchstim = find(varswitch);
%         if (any(btdstruct.btd.findField([gq_switch '_ratePredROG_var_low']) <= 0) || any(btdstruct.btd.findField([gq_switch '_ratePredROG_var_low_last'])<=0) )
%             btdstruct = ratefunGQVariance(btdstruct, opstruct);
%         end
           gq_switch = [opstruct.stim(switchstim).gqname '_var'];
    otherwise
        warning ('I don''t know how to handle multiple variances switching simultaneously');
        return;
end
for j = 1:length(opstruct.stim)
    btdstruct = addScaledVarianceTurnLin (btdstruct, opstruct, j, gq_switch);
end

opstruct.switchstim = switchstim;
btdstruct.varops = opstruct;

btdstruct = gatherVarianceFields(btdstruct);

end

function btdstruct = addScaledVarianceTurnLin (btdstruct,opstruct,j, gq_switch)

gqname = opstruct.stim(j).gqname;
kernelTime = opstruct.kernelTime;
kernelDt = opstruct.kernelDt;

kname = [gqname 'TurnLin'];
btdstruct.linname{j} = kname;
btdstruct.kernelTime = kernelTime;
btdstruct.kernelDt = kernelDt;


[convkernel, btd] = btdstruct.btd.createBTAKernel('turn', 'start', gqname, kernelTime, kernelDt, 'newFieldName',kname, 'abbott', true);
stimInput = btd.behaviorTriggeredDataMatrix ('all', '', kname, 0);

existsAndDefault (gq_switch, []);
if (isempty(gq_switch) && opstruct.stim(j).switch)
    if (opstruct.stim(j).iswn)
        gq_switch = [opstruct.stim(j) '_var'];
    else
        gq_switch = opstruct.stim(j);
    end
end
if (~isempty(gq_switch))
    stimTon = btd.behaviorTriggeredDataMatrix ('all', '', [gq_switch '_ton'],0);
    stimToff = btd.behaviorTriggeredDataMatrix ('all', '', [gq_switch '_toff'],0);
    
    islow = stimToff < stimTon & stimToff > 0;
    inds = islow & stimToff > kernelTime; %guarantee only in low variance portion
    slow = std(stimInput(inds));
    ishigh = stimTon < stimToff & stimTon > 0;
    inds = ishigh & stimTon > kernelTime; %guarantee only in high variance portion
    shigh = std(stimInput(inds));
    s = min(slow, shigh);
    
else
    s = std(stimInput); %choose the lower variance epoch
end

convkernel = convkernel / s;
btdstruct.btd = btd.addConvolvedFields(gqname, kname, convkernel, kernelDt, 'scaleToSqr', true);
btdstruct.convkernels{j} = convkernel;

end

function btdstruct = gatherVarianceFields (btdstruct)

btd = btdstruct.btd;

switchstim = btdstruct.varops.switchstim;

aa = [btd.all]';
alleti = [aa.eti]';
tt = [btd.turn]';
turneti = [tt.start_eti]';
isrun = [aa.isrun]';


if (switchstim > 0)
    if (btdstruct.varops.stim(switchstim).iswn)
        gqswitch = [btdstruct.varops.stim(switchstim).gqname '_var'];
    else
        gqswitch = btdstruct.varops.stim(switchstim).gqname;
    end
    
    
    if(strcmpi(btdstruct.varops.stim(switchstim).ramptype, 'square'))
        gq_varh = strcat(gqswitch, '_high');
%         gq_varl = strcat(gqswitch, '_low');
    else
        gq_varh = strcat(gqswitch, '_rising');
%         gq_varl = strcat(gqswitch, '_falling');
    end
    
    ii = btd(1).findField(gq_varh);
    xd = btd(1).glt(ii).xData;
    yd = double(btd(1).glt(ii).yData);
    [freq, phase] = fitWave(xd, yd, 'square');
    period = 2*pi/freq;
    tshift = phase/freq;
    
    gq_ton = [gqswitch,'_ton'];
    gq_toff = [gqswitch,'_toff'];
    
    all_ton = btd.behaviorTriggeredDataMatrix('all', 'start', gq_ton, 0);
    all_toff = btd.behaviorTriggeredDataMatrix('all', 'start', gq_toff, 0);
    turn_ton = btd.behaviorTriggeredDataMatrix('turn', 'start', gq_ton, 0);
    turn_toff = btd.behaviorTriggeredDataMatrix('turn', 'start', gq_toff, 0);
else
    all_ton = [];
    all_toff = [];
    turn_ton = [];
    turn_toff = [];
end
%
%     if (isempty(btdstruct.varops.stim(switchstim).period))
%         btdstruct.varops.stim(switchstim).period = max(max(all_ton), max(all_toff));
%     end
%     period = btdstruct.varops.stim(switchstim).period;
mintime = min(alleti(all_ton > 0 & all_toff > 0));
btdstruct.varops.timeRange(1) = max(btdstruct.varops.timeRange(1), mintime);
turnvalid = turneti >= min(btdstruct.varops.timeRange) & turneti < max(btdstruct.varops.timeRange);
allvalid = alleti >= min(btdstruct.varops.timeRange) & alleti < max(btdstruct.varops.timeRange);

for whichStim = find([btdstruct.varops.stim.iswn])
    
    
    gq = btdstruct.varops.stim(whichStim).gqname;
    gq_conv = strcat(gq, 'TurnLin');
    [turndata,~,texpnum] = btd.behaviorTriggeredDataMatrix('turn', 'start', gq_conv, 0);
    [alldata,~,aexpnum] =  btd.behaviorTriggeredDataMatrix('all', 'start', gq_conv, 0);
    btdstruct.var(whichStim).gqname =  btdstruct.varops.stim(whichStim).gqname;
    btdstruct.var(whichStim).fullensemble.x_conv = alldata(allvalid);
    btdstruct.var(whichStim).fullensemble.eti = alleti (allvalid);
    btdstruct.var(whichStim).fullensemble.ton = all_ton (allvalid);
    btdstruct.var(whichStim).fullensemble.toff = all_toff (allvalid);
    btdstruct.var(whichStim).fullensemble.expnum = aexpnum (allvalid);
    
    btdstruct.var(whichStim).noturn.x_conv = alldata(isrun & allvalid);
    btdstruct.var(whichStim).noturn.eti = alleti (isrun & allvalid);
    btdstruct.var(whichStim).noturn.ton = all_ton (isrun & allvalid);
    btdstruct.var(whichStim).noturn.toff = all_toff (isrun & allvalid);
    btdstruct.var(whichStim).noturn.expnum = aexpnum  (isrun & allvalid);
    
    btdstruct.var(whichStim).turn.x_conv = turndata(turnvalid);
    btdstruct.var(whichStim).turn.eti = turneti (turnvalid);
    btdstruct.var(whichStim).turn.ton = turn_ton (turnvalid);
    btdstruct.var(whichStim).turn.toff = turn_toff (turnvalid);
    btdstruct.var(whichStim).turn.expnum = texpnum (turnvalid);
    btdstruct.var(whichStim).period = period;
    btdstruct.var(whichStim).tshift = tshift;
end



end




