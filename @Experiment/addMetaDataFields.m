function addMetaDataFields(expt, varargin)
%parses the metadata file (.mdat) for global quantities

if (length(expt) > 1)  
    for j = 1:length(expt)
        expt(j).addMetaDataFields(varargin{:});
    end
    return;
end
        
if (~isempty(varargin) && exist(varargin{1},'file'))
    expt.timfname = varargin{1};
end

% header = {'GAS_CH1_SP',  'co2filt_0', 'co2filt_1', 'PID_PPM', 'LED_ON', 'LUM_PWM_AMP'};
% timheader = {'GAS_CH1_SP_time', 'co2fast_0_time', 'co2fast_1_time', 'PID_PPM_time', 'LED_ON_time', 'LUM_PWM_AMP_time'};
% gname = {'gassp', 'co2ppm0', 'co2ppm1', 'vocppm', 'ledOn', 'ledAmp'};
% smoothTime = [0, 2, 2, 7, 0, 0];
% derivTime = [0.1, 1, 1, 1, 0, 0];
% isgas = [true, true, true, true, false, false];
header = {'GAS_CH1_SP',  'co2filt_0', 'co2filt_1', 'PID_PPM', 'LED_ON', 'LUM_PWM_AMP', 'Temp_C_WS'};
timheader = {'GAS_CH1_SP_time', 'co2fast_0_time', 'co2fast_1_time', 'PID_PPM_time', 'LED_ON_time', 'LUM_PWM_AMP_time','Temp_C_WS_time'};
gname = {'gassp', 'co2ppm0', 'co2ppm1', 'vocppm', 'ledOn', 'ledAmp','temperature'};
smoothTime = [0, 2, 2, 7, 0, 0,1];
derivTime = [0.1, 1, 1, 1, 0, 0,1];
isgas = [true, true, true, true, false, false,false];

gasvel = 2.5; %changed to 2.5 for 4L/min flow in new flow chamber 6/12/2016
gasdir = [-1;0]; %new flow chamber, flow is from the right 
gasorigin_cell = {[27;0],[27;0],[27;0],[27;0],[],[],[]}; %changed to +27 cm for new flow chamber 6/12/2016
fnum = 1:length(expt.elapsedTime);
et = expt.elapsedTime;

if (isempty(expt.metadata))
    expt.metadata = importdata2(fixFileNameWin(expt.timfname));
end

datastruct = expt.metadata;

%{
try
    datastruct = importdata(timfname);
catch me
    disp(me);
    return;
end
%}

for j = 1:length(header)
    col = find(strcmpi(datastruct.colheaders, header{j}), 1);
    if (~isempty(col))
        try
            ydata = datastruct.data(:,col);
            col = find(strcmpi(datastruct.colheaders, timheader{j}), 1);
            if (~isempty(col))
                timdata = datastruct.data(:,col)/1000;
            else
                timdata = et;
            end
            inds = find(isfinite(ydata)&isfinite(timdata));
            ydat = ydata(inds);
            timdata = timdata(inds);
            [~,I] = min(timdata);
            ydat = ydat(I:end); %get rid of any early bad time values before reset
            timdata = timdata(I:end);
            [timdata,I] = unique(timdata);
            ydat = ydat(I);
            p = polyfit(timdata, et(inds(I)),1);
            et2 = p(1)*timdata + p(2);
            eti = linspace(min(et2), max(et2), (max(et2)-min(et2))/expt.dr.interpTime);
            xdat = eti;
            ydat = interp1(et2, ydat, eti, 'linear');
            inds = isfinite(xdat) & isfinite(ydat);
            xdat = xdat(inds); ydat = ydat(inds);
            if (smoothTime(j) > 0)
                ydat = lowpass1D(ydat, smoothTime(j)/expt.dr.interpTime, 'padType', 'linear');
            end
            gq = GlobalQuantity();
            gq.yData = ydat;
            gq.fieldname = gname{j};
            if (isgas(j))
                gasorigin = gasorigin_cell{j};
                gq.xField = {'eti', 'sloc'};
                xData.et = xdat;
                xData.origin = gasorigin;
                xData.flowspeed = gasvel;
                xData.flowdir = gasdir;
                gq.xData = xData;
                gq.derivationMethod = @GlobalQuantity.timeVaryingGasDerivation;
            else
                gq.xField = 'eti';
                gq.xData = [min(expt.elapsedTime)-1 xdat max(expt.elapsedTime)+1];
                gq.yData = gq.yData([1 1:end end]);
            end
        
            
        
            expt.addGlobalQuantity(gq);
            if (derivTime(j) > 0)
                gq.xData = [];
                gq.yData = [];
                gq.xField = gname{j};
                gq.fieldname = ['d' gname{j}];
                gq.derivationMethod = @(xin, xData, yData) deriv(xin, derivTime(j)/expt.dr.interpTime);
                expt.addGlobalQuantity(gq);
            end
        catch me
            disp (['problem adding metadata field ' gname{j}]);
            disp (me.getReport);
        end
            
    end
end

ind = find(strcmpi('LAST_STROBE_MS', datastruct.colheaders));
if (~isempty(ind))
    strobeind = find(isfinite(datastruct.data(:,ind)));
    strobeind = strobeind(strobeind > 0);
    xdata = et;
    ydata = false(size(xdata));
    ydata(strobeind) = true;
    
    gq = GlobalQuantity();
    gq.yData = ydata;
    gq.fieldname = 'strobe';
    gq.xField = 'eti';
    gq.xData = xdata;
         
    expt.addGlobalQuantity(gq);
    gq.fieldname = 'strobe_high';
    expt.addGlobalQuantity(gq);
    gq.fieldname = 'strobe_low';
    gq.yData = ~ydata;
    expt.addGlobalQuantity(gq);
    
    
    strobetime = et(strobeind);
    ydata = -100*ones(size(xdata));
    for j = 1:length(et)
        if (any(strobetime <= et(j)))
            ydata(j) = et(j) - max(strobetime(strobetime <= et(j)));
        end
    end
    xdata = et;
    gq.yData = ydata;
    gq.fieldname = 'strobe_ton';
    gq.xField = 'eti';
    gq.xData = xdata;
    expt.addGlobalQuantity(gq);
    gq.fieldname = 'strobe_toff';
    expt.addGlobalQuantity(gq);
end

if (isfield(expt.metadata, 'ledNBitsOut') || isfield(expt.metadata, 'ledNBytesOut'))
    [p,f] = fileparts(expt.fname);
    if (isfield(expt.metadata, 'ledNBitsOut'))
        nb = expt.metadata.ledNBitsOut;
        fp = expt.metadata.ledFP;
        ii = fp > 0 & isfinite(nb ./ fp); 
        if (median(nb(ii)./fp(ii)) < 4) %byte output, not bit output
            expt.metadata.ledNBytesOut = expt.metadata.ledNBitsOut;
            expt.metadata = rmfield(expt.metadata, 'ledNBitsOut');
        end 
    end
    
    try
        createLedTablesFromBinFiles(expt);
%         fname = fullfile (p, [f ' led values.bin']);
%         if (exist(fname, 'file'))
%             [glt,gq] = GlobalLookupTable.createLedTableFromBitFile(expt, fname, true); %true means add global quantities to experiment
%         else
%             fname = fullfile (p, [f ' sup data dir'], [f ' led values.bin']);
%             if (exist(fname, 'file'))
%                 [glt,gq] = GlobalLookupTable.createLedTableFromBitFile(expt, fname, true);
%             else
%                 glt = [];
%             end
%         end
    catch me
        glt = [];
        disp (me.getReport());
    end
    %gqs already added
    %{
    if (~isempty(glt))
        expt.addGlobalQuantity(glt);
        expt.addGlobalQuantity(gq);
    end
    %} 
end
    
end

function createLedTablesFromBinFiles(expt)
    [p,f] = fileparts(expt.fname);
    fend = {' led values.bin', ' led1 values.bin', ' led2 values.bin', '@led1.bin', '@led2.bin'};
    gqname = {'ledVal', 'led1Val','led2Val', 'led1Val','led2Val'};
    found = false(size(gqname));
    for j = 1:length(fend)
        fname = fullfile (p, [f fend{j}]);

        if (exist(fixFileNameWin(fname), 'file'))
            [glt,gq] = GlobalLookupTable.createLedTableFromBitFile(expt, fname, false); %true means add global quantities to experiment
        else
            sdd = Experiment.findSuppDataDir(expt.fname);
            fname = fullfile (sdd, [f fend{j}]);
            if (~isempty(sdd) && exist(fixFileNameWin(fname), 'file'))
                [glt,gq] = GlobalLookupTable.createLedTableFromBitFile(expt, fname, false);
            else
                glt = [];
            end
        end
        if (~isempty(glt))
            found(j) = true;
            glt(1).fieldname = gqname{j};
            glt(2).fieldname = [gqname{j} 'Deriv'];
            glt(3).fieldname = [gqname{j} 'Diff'];
            gq(1).fieldname = gqname{j};
            gq(2).fieldname = [gqname{j} 'Deriv'];
            gq(3).fieldname = [gqname{j} 'Diff'];
            expt.addGlobalQuantity(glt);
            expt.addGlobalQuantity(gq);
        end
    end
    if (~any(found))
        for j = 1:length(fend)
            %relax rules for matching; see if there is anything with the right
            %ending in sup data dir
            sdd = Experiment.findSuppDataDir(expt.fname);
            if (isempty(sdd))
                continue;
            end
            dd = dir(fullfile(sdd, ['*' fend{j}]));
            if (isempty(dd))
                continue;
            end
            if (length(dd) == 1)
                fname = fullfile(sdd, dd.name);
            else
                names = {dd.name};
                for n = 1:length(f)
                    match = strncmpi(f, names, n);
                    if (nnz(match) == 1)
                        fname = fullfile(sdd, dd(match).name);
                        break;
                    end
                    if (nnz(match) == 0)
                        warning ('expt:amd', ['multiple equally good/poor choices for '  fend{j} ' for: ' f]);
                        continue;
                    end
                end
            end
            [glt,gq] = GlobalLookupTable.createLedTableFromBitFile(expt, fname, false);
            if (~isempty(glt))
                found(j) = true;
                glt(1).fieldname = gqname{j};
                glt(2).fieldname = [gqname{j} 'Deriv'];
                glt(3).fieldname = [gqname{j} 'Diff'];
                gq(1).fieldname = gqname{j};
                gq(2).fieldname = [gqname{j} 'Deriv'];
                gq(3).fieldname = [gqname{j} 'Diff'];
                expt.addGlobalQuantity(glt);
                expt.addGlobalQuantity(gq);
            end
                        
        end
    end
    if (~any(found))
        warning('expt:amd', [' can not find led value file for '  f]);
    end
    
        

end        
        
