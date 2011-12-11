function gqs = timeOnOffGQs (gq, ramptype, varargin)
%function gqs = timeOnOffGQs (gq, ramptype, varargin)
%
%creates time on/off fields based on a global quantity
%the time on is the time since the quantity went high (or started rising
%for a ramp)
%
%ramptype is 'square', 'triangle', 'sine', 'exponential', 

if (isstruct(gq.xData))
    et = gq.xData.et;
else
    et = gq.xData;
end
ydat = gq.yData;
gqs = repmat(gq, [1 4]);

switch (lower(ramptype))
    case 'square'
        [idx, c] = kmeans(ydat', 2, 'start', [min(ydat);max(ydat)]);
        if (c(1) > c(2))
            high = idx' == 1;
            low = idx' == 2;
        else
            high = idx' == 2;
            low = idx' == 1;
        end
        
        high = medfilt2(high, [1 9], 'symmetric');
        low = medfilt2(low, [1 9], 'symmetric');
        gqs(1).fieldname = [gq.fieldname '_high'];
        gqs(2).fieldname = [gq.fieldname '_low'];        
        gqs(1).yData = high;
        gqs(2).yData = low;
        starthigh = find(diff(high) > 0);
        startlow = find(diff(low) > 0);
        dv = deriv(ydat, diff(starthigh(1:2)/20));
        valid = true(size(starthigh));        
        for j = 1:length(starthigh)
            ind = find(dv(1:starthigh(j)) < 0.1*dv(starthigh(j)), 1, 'last');
            if (~isempty(ind))
                starthigh(j) = ind +find(ydat(ind:end) > ydat(ind), 1, 'first');
            else
                valid(j) = false;
            end
        end
        starthigh = starthigh(valid);
        valid = true(size(startlow));        
       % figure(1);
      %  plot (1:length(dv), dv, 'b-', startlow, dv(startlow), 'g.', starthigh, dv(starthigh), 'r.');
        for j = 1:length(startlow)
            ind = find(dv(1:startlow(j)) > 0.1*dv(startlow(j)), 1, 'last');
            if (~isempty(ind))
                startlow(j) = ind +find(ydat(ind:end) < ydat(ind), 1, 'first');
            else
                valid(j) = false;
                
            end
        end
        startlow = startlow(valid);
     %   size(startlow)
    case {'triangle','exponential'}
        %{
        dv = deriv(ydat, 20);
        falling = medfilt2(dv < 0, [1 9], 'symmetric');
        rising = medfilt2(dv > 0, [1 9], 'symmetric');
        gqs(1).fieldname = [gq.fieldname '_rising'];
        gqs(2).fieldname = [gq.fieldname '_falling'];        
        gqs(1).yData = rising;
        gqs(2).yData = falling;
        %hx = linspace(min(dv), max(dv), 100);
        %figure(1); plot (hx, hist(dv(falling), hx), 'b-', hx, hist(dv(rising), hx), 'r-');
        %figure(2); plot (et(falling), ydat(falling), 'b.', et(rising), ydat(rising), 'r.');
        starthigh = find(diff(rising) > 0);
        startlow = find(diff(falling) > 0);
        %}
        starthigh = findPeriodicMaxima(-ydat);
        startlow = findPeriodicMaxima(ydat);
        
        falling = zeros(size(ydat));
        falling(startlow) = 1;
        falling(starthigh(starthigh > startlow(1))) = -1;
        falling = cumsum(falling);
        
        rising = zeros(size(ydat));
        rising(starthigh) = 1;
        rising(startlow(startlow > starthigh(1))) = -1;
        rising = cumsum(rising);
        
        gqs(1).fieldname = [gq.fieldname '_rising'];
        gqs(2).fieldname = [gq.fieldname '_falling'];        
        gqs(1).yData = logical(rising);
        gqs(2).yData = logical(falling);
        
    otherwise
        disp ('not implemented yet');
        return;
end

etstart = et(find(diff(et) < 10*median(diff(et)), 1, 'first'));
etend = et(find(diff(et) < 10*median(diff(et)), 1, 'last') + 1);

etf = [et(et < etstart) (etstart:(median(diff(et))/100):etend) et(et>etend)];
%length(etf)
%median(diff(etf))

ton = -100 * ones(size(etf));
toff = ton;
        
hightime = et(starthigh);
lowtime = et(startlow);

for j = 1:length(hightime)
    ton(etf >= hightime(j)) = etf(etf >= hightime(j)) - hightime(j);
end
for j = 1:length(lowtime)
    toff(etf >= lowtime(j)) = etf(etf >= lowtime(j)) - lowtime(j);
end
% for j = 1:length(startlow)
%     toff(startlow(j):end) = et(startlow(j):end) - et(startlow(j));
% end

xData = gqs(3).xData;
if (isstruct(xData))
    xData.et = etf;
else
    xData = etf;
end
gqs(3).xData = xData;
gqs(4).xData = xData;

if(strcmpi(func2str(gqs(3).derivationMethod), 'GlobalQuantity.oneDinterpolation'))
    gqs(3).derivationMethod = @GlobalQuantity.oneDinterpolationNearest;
    gqs(4).derivationMethod = gqs(3).derivationMethod;
end
if(strcmpi(func2str(gqs(3).derivationMethod), 'GlobalQuantity.timeVaryingGasDerivation'))
    gqs(3).derivationMethod = @GlobalQuantity.timeVaryingGasDerivationNearest;
    gqs(4).derivationMethod = gqs(3).derivationMethod;
end

gqs(3).fieldname = [gq.fieldname '_ton'];
gqs(4).fieldname = [gq.fieldname '_toff'];
gqs(3).yData = ton;
gqs(4).yData = toff;

