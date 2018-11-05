function gqs = timeOnOffGQsVariance (gq, ramptype, varargin)
%function gqs = timeOnOffGQsVariance (gq, ramptype, varargin)
%
%creates time on/off fields based on the magnitude of the variance of a 
%particular gq 
%calculates variance, discarding outliers then calls timeOnOffGQsVariance
%see timeOnOffGQs for more information


ydat = gq.yData;

%discard outliers/transients and find variance
yds = (ydat-median(ydat)).^2;
valid = (yds < 10*percentile(yds,.99));
yds = (ydat - mean(ydat(valid))).^2;
yds(~valid) = 0; %variance

gq.yData = yds;
gq.fieldname = [gq.fieldname '_var'];

gqs = gq.timeOnOffGQs(ramptype, varargin{:});

return
%old code - eliminate when sure new code works
sigma = 0.2 / median(diff(et)); % average variance over 0.2 seconds
ydat = medfilt1(lowpass1D(yds, sigma), ceil(sigma)); %median filter over 0.2 second



gqs = repmat(gq, [1 6]);

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
        
        high = ~notany&medfilt2(high, [1 ceil(1/median(diff(et)))], 'symmetric');
        low = ~notany&medfilt2(low, [1 ceil(1/median(diff(et)))], 'symmetric');
        gqs(1).fieldname = [gq.fieldname '_var_high'];
        gqs(2).fieldname = [gq.fieldname '_var_low'];        
        gqs(1).yData = high;
        gqs(2).yData = low;
        %this marks the point before the transition
        starthigh = find(diff(high) > 0 &  ~(diff(notany) < 0));
        startlow = find(diff(low) > 0 & ~(diff(notany) < 0));
%         oldstarthigh = starthigh;
%         dv = deriv(ydat, median(diff(starthigh))/20); %this is a disgusting kludge!  fix it.  why should the derivative be related to the period?
%         valid = true(size(starthigh));        
%         for j = 1:length(starthigh)
%             %find the last point before the derivative starts increasing
%             ind = find(dv(1:starthigh(j)) < 0.1*dv(starthigh(j)), 1, 'last');
%             if (~isempty(ind))
%                 %find the first point with a value above the value at that
%                 %derivative point (shouldn't it always be the next point if
%                 %derivative is positive?)
%                 starthigh(j) = ind +find(ydat(ind:end) > ydat(ind), 1, 'first') -1;
%             else
%                 valid(j) = false;
%             end
%         end
%         starthigh = starthigh(valid);
% %         figure(1);
% %         plot (gqs(1).xData, high, gqs(1).xData(oldstarthigh),high(oldstarthigh), 'r.', gqs(1).xData(starthigh), high(starthigh), 'g.');
% %         figure(2);
% %         plot (gqs(1).xData, dv, gqs(1).xData(oldstarthigh),dv(oldstarthigh), 'r.', gqs(1).xData(starthigh), dv(starthigh), 'g.');
% %         
%         valid = true(size(startlow));        
%        % figure(1);
%       %  plot (1:length(dv), dv, 'b-', startlow, dv(startlow), 'g.', starthigh, dv(starthigh), 'r.');
%         for j = 1:length(startlow)
%             ind = find(dv(1:startlow(j)) > 0.1*dv(startlow(j)), 1, 'last');
%             if (~isempty(ind))
%                 startlow(j) = ind +find(ydat(ind:end) < ydat(ind), 1, 'first') - 1;
%             else
%                 valid(j) = false;
%                 
%             end
%         end
%         startlow = startlow(valid);
     %   size(startlow)
        if (false)
            plot (et(high), ydat(high), 'r.', et(low), ydat(low), 'b.');
            pause;
        end
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

gqs(3).fieldname = [gq.fieldname '_var_ton'];
gqs(4).fieldname = [gq.fieldname '_var_toff'];
gqs(3).yData = ton;
gqs(4).yData = toff;

gqs(5).fieldname = [gq.fieldname '_var_cyclenum_on'];
gqs(6).fieldname = [gq.fieldname '_var_cyclenum_off'];
gqs(5).yData = whichon;
gqs(6).yData = whichoff;





