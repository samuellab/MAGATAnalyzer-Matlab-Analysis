function btdstruct = prepVarianceSwitchingAnalysis_Gepner(btdstruct, opstruct)
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
        if(strcmpi(opstruct.stim(j).gqname, 'led2m1ValDiff')) 
            if(strcmpi(opstruct.stim(1).gqname, 'led1ValDiff'))
                stims{1} = 'led1Val';
                stimdiffs{1} = 'led1ValDiff';
            else
                stims{1} = 'ledVal';
                stimdiffs{1} = 'ledValDiff';
            end
            if(strcmpi(opstruct.stim(2).gqname, 'led2ValDiff'))
                stims{2} = 'led2Val';
                stimdiffs{2} = 'led2ValDiff';
            else
                stims{2} = 'ledVal';
                stimdiffs{2} = 'ledValDiff';
            end
            OdorPlusLight = @(yd) 1/sqrt(2) * (yd{1} - yd{2});
            btdstruct.btd = btdstruct.btd.addOperationFields({stims{2}, stims{1}}, 'led2m1Val', OdorPlusLight);
            btdstruct.btd = btdstruct.btd.addOperationFields({stimdiffs{2}, stimdiffs{1}}, opstruct.stim(j).gqname, OdorPlusLight);
        end
        if (any(btdstruct.btd.findField([opstruct.stim(j).gqname '_var_ton']) <= 0))
            btdstruct.btd = btdstruct.btd.addVarianceGQs(opstruct.stim(j).gqname, opstruct.stim(j).ramptype, 'fixedPeriod', opstruct.stim(j).period);
        end
    end
    
end

for j = 1:length(opstruct.stim)
    %     if (opstruct.stim(j).iswn && ~strcmpi(opstruct.stim(j).gqname, 'led2m1ValDiff') && (opstruct.redogq || any(btdstruct.btd.findField([opstruct.stim(j).gqname 'TurnLin']) <= 0)))
    if ( ~strcmpi(opstruct.stim(j).gqname, 'led2m1ValDiff') && (opstruct.redogq || any(btdstruct.btd.findField([opstruct.stim(j).gqname 'TurnLin']) <= 0)))
        btdstruct = addScaledVarianceTurnLin (btdstruct, opstruct, j);
    end
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

opstruct.switchstim = switchstim;
btdstruct.varops = opstruct;

btdstruct = gatherVarianceFields(btdstruct, opstruct);

end

function [btdstruct, var] = gatherVarianceFields (btdstruct, opstruct)

btd = btdstruct.btd;

switchstim = btdstruct.varops.switchstim;

aa = [btd.all]';
alleti = [aa.eti]';
tt = [btd.turn]';
turneti = [tt.start_eti]';
isrun = [aa.isrun]';

acchs = [btd.acchs]';
rejhs = [btd.rejhs]';
acchseti = [acchs.start_eti]';
rejhseti = [rejhs.start_eti]';

accnhs = [acchs.hsnum]';
rejnhs = [rejhs.hsnum]';
firstacchs = accnhs == 1;
firstrejhs = rejnhs == 1;

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
    
    firsths.name = 'hsnum';
    firsths.validop = @(x) x == 1;

    all_ton = btd.behaviorTriggeredDataMatrix('all', 'start', gq_ton, 0);
    all_toff = btd.behaviorTriggeredDataMatrix('all', 'start', gq_toff, 0);
    turn_ton = btd.behaviorTriggeredDataMatrix('turn', 'start', gq_ton, 0);
    turn_toff = btd.behaviorTriggeredDataMatrix('turn', 'start', gq_toff, 0);
    acchs_ton = btd.behaviorTriggeredDataMatrix('acchs', 'start', gq_ton, 0);
    rejhs_ton = btd.behaviorTriggeredDataMatrix('rejhs', 'start', gq_ton, 0);
    acchs_toff = btd.behaviorTriggeredDataMatrix('acchs', 'start', gq_toff, 0);
    rejhs_toff = btd.behaviorTriggeredDataMatrix('rejhs', 'start', gq_toff, 0);
    
else
    all_ton = [];
    all_toff = [];
    turn_ton = [];
    turn_toff = [];
end

if (switchstim>0 && isempty(btdstruct.varops.stim(switchstim).period))
    btdstruct.varops.stim(switchstim).period = max(max(all_ton), max(all_toff));
    period = btdstruct.varops.stim(switchstim).period;
end
 

mintime = min(alleti(all_ton > 0 & all_toff > 0));
if(~isempty(mintime))
    btdstruct.varops.timeRange(1) = min(btdstruct.varops.timeRange(1), mintime);
end
    
turnvalid = turneti >= min(btdstruct.varops.timeRange) & turneti < max(btdstruct.varops.timeRange);
allvalid = alleti >= min(btdstruct.varops.timeRange) & alleti < max(btdstruct.varops.timeRange);
acchsvalid = acchseti >= min(btdstruct.varops.timeRange) & acchseti < max(btdstruct.varops.timeRange);
rejhsvalid = rejhseti >= min(btdstruct.varops.timeRange) & rejhseti < max(btdstruct.varops.timeRange);

for whichStim = 1:length(btdstruct.varops.stim)
    
    if(strcmpi(btdstruct.varops.stim(whichStim).gqname, 'led2m1ValDiff'))
        continue
    end
    
    gq = btdstruct.varops.stim(whichStim).gqname;
    gq_conv = strcat(gq, 'TurnLin');
    if(any(btdstruct.btd.findField(gq_conv) <= 0)) %for the case where gq is led1Val or led2Val
        [turndata,~,texpnum] = btd.behaviorTriggeredDataMatrix('turn', 'start', gq, 0);
        [alldata,~,aexpnum] =  btd.behaviorTriggeredDataMatrix('all', 'start', gq, 0);
    else
        [turndata,~,texpnum] = btd.behaviorTriggeredDataMatrix('turn', 'start', gq_conv, 0);
        [alldata,~,aexpnum] =  btd.behaviorTriggeredDataMatrix('all', 'start', gq_conv, 0);
    end
    var(whichStim).gqname =  btdstruct.varops.stim(whichStim).gqname;
    var(whichStim).fullensemble.x_conv = alldata(allvalid);
    var(whichStim).fullensemble.eti = alleti (allvalid);

    var(whichStim).fullensemble.expnum = aexpnum (allvalid);  
    var(whichStim).noturn.x_conv = alldata(isrun & allvalid);
    var(whichStim).noturn.eti = alleti (isrun & allvalid);

    var(whichStim).noturn.expnum = aexpnum  (isrun & allvalid);
    var(whichStim).turn.x_conv = turndata(turnvalid);
    var(whichStim).turn.eti = turneti (turnvalid);
    var(whichStim).turn.expnum = texpnum (turnvalid);

    var(whichStim).acchs.eti = acchseti(acchsvalid & firstacchs);
    var(whichStim).rejhs.eti = rejhseti(rejhsvalid & firstrejhs);
    
    if(switchstim > 0)
        var(whichStim).fullensemble.ton = all_ton (allvalid);
        var(whichStim).fullensemble.toff = all_toff (allvalid);
        var(whichStim).noturn.ton = all_ton (isrun & allvalid);
        var(whichStim).noturn.toff = all_toff (isrun & allvalid);
        var(whichStim).turn.ton = turn_ton (turnvalid);
        var(whichStim).turn.toff = turn_toff (turnvalid);
        
        var(whichStim).acchs.ton = acchs_ton(acchsvalid & firstacchs);
        var(whichStim).acchs.toff = acchs_toff(acchsvalid & firstacchs);
        var(whichStim).rejhs.ton = rejhs_ton(rejhsvalid & firstrejhs);
        var(whichStim).rejhs.toff = rejhs_toff(rejhsvalid & firstrejhs);
        
        var(whichStim).period = period;
        var(whichStim).tshift = tshift;
    end
    
    btdstruct.var = var;
end

if(switchstim>0)
    if(strcmpi(btdstruct.varops.stim(switchstim).gqname, 'led2m1ValDiff'))
        kname_uv = {'linU', 'linV'};
        
        stimnames = {opstruct.stim(1).gqname, opstruct.stim(2).gqname};
        kname1 = [stimnames{1} 'TurnLin'];
        kname2 = [stimnames{2} 'TurnLin'];

        if(isfield(opstruct, 'thetaUV') && ~isfield(opstruct.thetaUV, 'low') )
            
            theta_all = opstruct.thetaUV * ones(size(all_ton));
            theta_turn = opstruct.thetaUV * ones(size(turn_ton));
        
        elseif(isfield(btdstruct, 'theta_all'))
            theta_all = btdstruct.theta_all;
        elseif(isfield(opstruct, 'thetaUV') && isfield(opstruct.thetaUV, 'low') )
            
            if(strcmpi(opstruct.stim(switchstim).ramptype, 'square'))
                theta_all = zeros(size(all_ton));
                theta_all(all_ton<all_toff) = opstruct.thetaUV.high;
                theta_all(all_toff<all_ton) = opstruct.thetaUV.low;
                theta_all(theta_all==0) = pi/4;
                
                theta_turn = zeros(size(turn_ton));
                theta_turn(turn_ton<turn_toff) = opstruct.thetaUV.high;
                theta_turn(turn_toff<turn_ton) = opstruct.thetaUV.low;
                theta_turn(theta_turn==0) = pi/4;
            else
                
                
                theta_all = zeros(size(all_ton));
                theta_all(all_ton>(period/2-opstruct.adaptationTime) & all_ton<(period/2+opstruct.adaptationTime)) = opstruct.thetaUV.high;
                theta_all(all_toff>(period/2-opstruct.adaptationTime) & all_toff<(period/2+opstruct.adaptationTime)) = opstruct.thetaUV.low;
                if(isfield(opstruct.thetaUV, 'up'))
                    theta_all(all_ton>opstruct.adaptationTime) = opstruct.thetaUV.up;
                    theta_all(all_toff>opstruct.adaptationTime) = opstruct.thetaUV.down;
                end
                theta_all(theta_all==0) = pi/4;
                
                theta_turn = zeros(size(turn_ton));
                theta_turn(turn_ton>(period/2-opstruct.adaptationTime) & turn_ton<(period/2+opstruct.adaptationTime)) = opstruct.thetaUV.high;
                theta_turn(turn_toff>(period/2-opstruct.adaptationTime) & turn_toff<(period/2+opstruct.adaptationTime)) = opstruct.thetaUV.low;
                if(isfield(opstruct.thetaUV, 'up'))
                    theta_turn(turn_ton>opstruct.adaptationTime) = opstruct.thetaUV.up;
                    theta_turn(turn_toff>opstruct.adaptationTime) = opstruct.thetaUV.down;
                end
                theta_turn(theta_turn==0) = pi/4;
            end
            
        end
        
        alldata =  btd.behaviorTriggeredDataMatrix('all', 'start', {kname1, kname2}, 0);
        turndata =  btd.behaviorTriggeredDataMatrix('turn', 'start', {kname1, kname2}, 0);
        ad(:,1) = cos(theta_all).*alldata(:,1) + sin(theta_all).*alldata(:,2);
        td(:,1) = cos(theta_turn).*turndata(:,1) + sin(theta_turn).*turndata(:,2);
        ad(:,2) = cos(theta_all).*alldata(:,2) - sin(theta_all).*alldata(:,1);
        td(:,2) = cos(theta_turn).*turndata(:,2) - sin(theta_turn).*turndata(:,1);
        
        for whichStim=1:2
            
%             [td,~,texpnum] = btdstruct.btd.behaviorTriggeredDataMatrix('turn', 'start', kname_uv{whichStim}, 0);
%             [ad,~,aexpnum] = btdstruct.btd.behaviorTriggeredDataMatrix('all', 'start', kname_uv{whichStim}, 0);
            btdstruct.var_uv(whichStim).gqname = kname_uv{whichStim};

            btdstruct.var_uv(whichStim).fullensemble.theta = theta_all(allvalid);
            btdstruct.var_uv(whichStim).fullensemble.x_conv = ad(allvalid, whichStim);
            btdstruct.var_uv(whichStim).fullensemble.eti = alleti (allvalid);
            btdstruct.var_uv(whichStim).fullensemble.ton = all_ton (allvalid);
            btdstruct.var_uv(whichStim).fullensemble.toff = all_toff (allvalid);
            btdstruct.var_uv(whichStim).fullensemble.expnum = aexpnum (allvalid);
            
            btdstruct.var_uv(whichStim).noturn.theta = theta_all(isrun & allvalid);
            btdstruct.var_uv(whichStim).noturn.x_conv = ad(isrun & allvalid, whichStim);
            btdstruct.var_uv(whichStim).noturn.eti = alleti (isrun & allvalid);
            btdstruct.var_uv(whichStim).noturn.ton = all_ton (isrun & allvalid);
            btdstruct.var_uv(whichStim).noturn.toff = all_toff (isrun & allvalid);
            btdstruct.var_uv(whichStim).noturn.expnum = aexpnum  (isrun & allvalid);
            
            btdstruct.var_uv(whichStim).turn.theta = theta_turn(turnvalid);
            btdstruct.var_uv(whichStim).turn.x_conv = td(turnvalid, whichStim);
            btdstruct.var_uv(whichStim).turn.eti = turneti (turnvalid);
            btdstruct.var_uv(whichStim).turn.ton = turn_ton (turnvalid);
            btdstruct.var_uv(whichStim).turn.toff = turn_toff (turnvalid);
            btdstruct.var_uv(whichStim).turn.expnum = texpnum (turnvalid);
            btdstruct.var_uv(whichStim).period = period;
            btdstruct.var_uv(whichStim).tshift = tshift;
            
        end
    end
end

% for whichStim = find([btdstruct.varops.stim.iswn])




end

function btdstruct = addScaledVarianceTurnLin (btdstruct,opstruct,j)

switchstim = find([opstruct.stim.switch], 1, 'first');
if(~isempty(switchstim))
    switchname = opstruct.stim(switchstim).gqname;
end
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

 if( (nnz([opstruct.stim.iswn]) == 1 && opstruct.stim(j).switch == 1 && strcmpi(opstruct.stim(j).ramptype, 'square')) || nnz([opstruct.stim.switch])==0 || (nnz([opstruct.stim.iswn]) == 2 && ~isfield(opstruct, 'convkernels')) )
%if( (nnz([opstruct.stim.iswn]) == 1 && opstruct.stim(j).switch == 1 && strcmpi(opstruct.stim(j).ramptype, 'square')) || (nnz([opstruct.stim.iswn]) == 2 && ~isfield(opstruct, 'convkernels')) )

    if(~isempty(switchstim))
        stimTon = btd.behaviorTriggeredDataMatrix ('all', '', [switchname '_var_ton'],0);
        stimToff = btd.behaviorTriggeredDataMatrix ('all', '', [switchname '_var_toff'],0);
        turnTon = btd.behaviorTriggeredDataMatrix ('turn', 'start', [switchname '_var_ton'],0);
        turnToff = btd.behaviorTriggeredDataMatrix ('turn', 'start', [switchname '_var_toff'],0);
        
        islow = turnToff < turnTon & turnToff > 0;
        ishigh = turnToff > turnTon & turnToff > 0;
        islow_all = stimToff < stimTon & stimToff > 0;
        ishigh_all = stimToff > stimTon & stimToff > 0;
    end
    %using only one kernel
    
    btdstruct.linname{j} = kname;
    [convkernel, btd] = btdstruct.btd.createBTAKernel('turn', 'start', gqname, kernelTime, kernelDt, 'newFieldName',kname, 'abbott', false);
    stimInput = btd.behaviorTriggeredDataMatrix ('all', '', kname, 0);
    turnInput = btd.behaviorTriggeredDataMatrix('turn', 'start', gqname, taxis);
    
    if(opstruct.stim(j).switch)
        inds_low = islow_all & stimToff > kernelTime; %guarantee only in low variance portion
        inds_high = ishigh_all & stimTon>kernelTime;
        sl = std(stimInput(inds_low));
        sh = std(stimInput(inds_high));
        s = sl;
        btdstruct.kernelScaling(j) = s;
    else
        s = std(stimInput) / opstruct.stim(j).var;
%         s = std(stimInput);
    end
    
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

    convkernel = opstruct.convkernels{j};

    
end

btdstruct.btd = btdstruct.btd.addConvolvedFields(gqname, kname, convkernel, kernelDt, 'scaleToSqr', true);

btdstruct.convkernels{j} = convkernel;


end



