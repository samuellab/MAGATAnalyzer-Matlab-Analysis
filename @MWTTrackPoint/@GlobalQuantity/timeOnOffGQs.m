function gqs = timeOnOffGQs (gq, ramptype, varargin)
%function gqs = timeOnOffGQs (gq, ramptype, varargin)
%
%creates time on/off fields based on a global quantity
%the time on is the time since the quantity went high (or started rising
%for a ramp)
%
%ramptype is 'square', 'triangle', 'sine', 'exponential',
%varargin - 'oldstyle', true uses method in place before 3/24/2016 update
%         - 'fixedPeriod', t - if the period is known, input it here. (not
%                               used by oldstyle)


if (isstruct(gq.xData))
    et = gq.xData.et;
else
    et = gq.xData;
end
ydat = gq.yData;
gqs = repmat(gq, [1 6]);

oldstyle = false;
varagin = assignApplicable(varargin);

if(~oldstyle)
    inds = ceil(.01*length(et)):floor(.99*length(et)); %discard outer 1% for fitting purposes
    inds = inds(inds-min(et) > 10); %and discard the first 10 seconds
    [freq, phase] = fitWave(et(inds), ydat(inds), ramptype, varargin{:});
    period = 2*pi/freq;
    ton = mysawtooth(freq*et + phase, 1, 0, period);
    toff = mysawtooth(freq*et + phase + pi, 1, 0, period);
    high = logical (sign(sin(freq*et + phase)) > 0);
    low = logical (sign(sin(freq*et + phase)) < 0);
    %
    %     rising = diff(yfit) > 0; rising(end+1) = rising(end);
    %     falling = diff(yfit) < 0; falling(end+1) = falling(end);
    onstart = [diff(ton) 1] < 0;
    offstart = [diff(toff) 1] < 0;
    cycleon = cumsum(onstart);
    cycleoff = cumsum(offstart);
    
    ton(cycleon < 1) = -100;
    toff(cycleoff < 1) = -100;
    
    switch (lower(ramptype))
        case 'square'
            gqs(1).fieldname = [gq.fieldname '_high'];
            gqs(2).fieldname = [gq.fieldname '_low'];
        otherwise
            gqs(1).fieldname = [gq.fieldname '_rising'];
            gqs(2).fieldname = [gq.fieldname '_falling'];
    end
    gqs(1).yData = high;
    gqs(2).yData = low;
    
    gqs(3).fieldname = [gq.fieldname '_ton'];
    gqs(4).fieldname = [gq.fieldname '_toff'];
    gqs(3).yData = ton;
    gqs(4).yData = toff;
    
    gqs(5).fieldname = [gq.fieldname '_cyclenum_on'];
    gqs(6).fieldname = [gq.fieldname '_cyclenum_off'];
    gqs(5).yData = cycleon;
    gqs(6).yData = cycleoff;
else
    
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
            
            nclust = 2;
            
            while (nnz(high) > 5*nnz(low) || nnz(low) > 5*nnz(high))
                nclust = nclust + 1;
                %bad clustering if we end up with too many of one kind
                if (nclust > 5)
                    error ('couldn''t partition data into high and low');
                end
                yx = linspace(min(ydat), max(ydat), 100);
                h1 = hist(ydat, yx);
                [~,I] = sort(h1, 'descend');
                %yx(I(1:nclust))
                
                [idx, c] = kmeans(ydat', nclust, 'start', yx(I(1:nclust))','emptyaction', 'singleton');%linspace(min(ydat),max(ydat),nclust)'
                nc = zeros([1 nclust]);
                for j = 1:nclust
                    nc(j) = nnz(idx == j);
                end
                [~,I] = sort(nc, 'descend');
                if (c(I(1)) > c(I(2)))
                    high = idx' == I(1);
                    low = idx' == I(2);
                else
                    high = idx' == I(2);
                    low = idx' == I(1);
                end
                
                
                c = c([I(1) I(2)]);
            end
            if (false)
                disp(['cluster centers are ' num2str(c')]);
            end
            notany = (~high & ~low);
            notany = medfilt2(notany, [1 9] , 'symmetric');
            notany = imdilate(notany, ones([1 5]));
            
            high = ~notany&medfilt2(high, [1 9], 'symmetric');
            low = ~notany&medfilt2(low, [1 9], 'symmetric');
            gqs(1).fieldname = [gq.fieldname '_high'];
            gqs(2).fieldname = [gq.fieldname '_low'];
            gqs(1).yData = high;
            gqs(2).yData = low;
            %this marks the point before the transition
            starthigh = find(diff(high) > 0 &  ~(diff(notany) < 0));
            startlow = find(diff(low) > 0 & ~(diff(notany) < 0));
            oldstarthigh = starthigh;
            dv = deriv(ydat, median(diff(starthigh))/20); %this is a disgusting kludge!  fix it.  why should the derivative be related to the period?
            valid = true(size(starthigh));
            for j = 1:length(starthigh)
                %find the last point before the derivative starts increasing
                ind = find(dv(1:starthigh(j)) < 0.1*dv(starthigh(j)), 1, 'last');
                if (~isempty(ind))
                    %find the first point with a value above the value at that
                    %derivative point (shouldn't it always be the next point if
                    %derivative is positive?)
                    starthigh(j) = ind +find(ydat(ind:end) > ydat(ind), 1, 'first') -1;
                else
                    valid(j) = false;
                end
            end
            starthigh = starthigh(valid);
            %         figure(1);
            %         plot (gqs(1).xData, high, gqs(1).xData(oldstarthigh),high(oldstarthigh), 'r.', gqs(1).xData(starthigh), high(starthigh), 'g.');
            %         figure(2);
            %         plot (gqs(1).xData, dv, gqs(1).xData(oldstarthigh),dv(oldstarthigh), 'r.', gqs(1).xData(starthigh), dv(starthigh), 'g.');
            %
            valid = true(size(startlow));
            % figure(1);
            %  plot (1:length(dv), dv, 'b-', startlow, dv(startlow), 'g.', starthigh, dv(starthigh), 'r.');
            for j = 1:length(startlow)
                ind = find(dv(1:startlow(j)) > 0.1*dv(startlow(j)), 1, 'last');
                if (~isempty(ind))
                    startlow(j) = ind +find(ydat(ind:end) < ydat(ind), 1, 'first') - 1;
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
    
    detf = max(median(diff(et))/100, 1E-3); %finest resolution we desire is 1 ms
    
    etf = [et(et < etstart) (etstart:detf:etend) et(et>etend)];
    %length(etf)
    %median(diff(etf))
    
    ton = -100 * ones(size(etf));
    toff = ton;
    whichon = -1 * ones(size(etf));
    whichoff = -1 * ones(size(etf));
    
    hightime = et(starthigh);
    lowtime = et(startlow);
    
    for j = 1:length(hightime)
        ton(etf >= hightime(j)) = etf(etf >= hightime(j)) - hightime(j);
        whichon(etf >= hightime(j)) = j;
    end
    for j = 1:length(lowtime)
        toff(etf >= lowtime(j)) = etf(etf >= lowtime(j)) - lowtime(j);
        whichoff(etf >= lowtime(j)) = j;
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
    for k = 3:6
        gqs(k).xData = xData;
    end
    % gqs(4).xData = xData;
    % gqs(5).xData = xData;
    % gqs(6).xData = xData;
    
    if(strcmpi(func2str(gqs(3).derivationMethod), 'GlobalQuantity.oneDinterpolation'))
        gqs(3).derivationMethod = @GlobalQuantity.oneDinterpolationNearest;
        for k = 4:6
            gqs(k).derivationMethod = gqs(3).derivationMethod;
        end
    end
    if(strcmpi(func2str(gqs(3).derivationMethod), 'GlobalQuantity.timeVaryingGasDerivation'))
        gqs(3).derivationMethod = @GlobalQuantity.timeVaryingGasDerivationNearest;
        for k = 4:6
            gqs(k).derivationMethod = gqs(3).derivationMethod;
        end
    end
    
    gqs(3).fieldname = [gq.fieldname '_ton'];
    gqs(4).fieldname = [gq.fieldname '_toff'];
    gqs(3).yData = ton;
    gqs(4).yData = toff;
    
    gqs(5).fieldname = [gq.fieldname '_cyclenum_on'];
    gqs(6).fieldname = [gq.fieldname '_cyclenum_off'];
    gqs(5).yData = whichon;
    gqs(6).yData = whichoff;
end