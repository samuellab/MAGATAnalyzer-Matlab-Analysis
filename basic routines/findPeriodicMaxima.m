function inds = findPeriodicMaxima(ydata, varargin)
%function inds = findPeriodicMaxima(ydata, varargin)
%
%finds maxima of the quasi-periodic function ydata, looking for robust
%maxima spaced somewhat evenly apart
debug = false;
varargin = assignApplicable(varargin);
ydata = ydata - mean(ydata);


%determine the period using autocorrelation
[xc,lags] = xcorr(ydata,ydata);
xc = xc(lags > 0);
lags = lags(lags > 0);

dxc = deriv(xc, 3);
firstmin = find(diff(sign(dxc)) > 0, 1, 'first');
allmax = find(diff(sign(dxc)) < 0);

firstmax = allmax(find(allmax > firstmin, 1, 'first'));

period = lags(firstmax);

%determine the phase

maxlen = period*floor(length(ydata) / period);
periodicmean = mean(reshape(ydata(1:maxlen), period, []),2);

[~,phase] = max(periodicmean);

%find low-passed local maxima within 0.1 period of the periodic maxima

yd = deriv(ydata, max(3,period/300));
zcr = find(diff(sign(yd)) < 0);
zcr = zcr(abs((zcr - phase)/period - round((zcr-phase)/period)) < 0.1);

if (debug)
    phase 
    period
    x = 1:length(yd);
    iv = (abs((x - phase)/period - round((x-phase)/period)) < 0.1);
   % plot (x, (x-phase)/period, x, round((x-phase)/period), x, abs((x - phase)/period - round((x-phase)/period)));
    plot (1:length(yd), yd, x(iv), yd(iv), 'g-', x(zcr), yd(zcr), 'r.');
end


valid = true(size(zcr));

for j = 1:length(zcr)
    if (nnz(abs(zcr - zcr(j)) < period/4) > 1)
        inds = find(abs(zcr - zcr(j)) < period/4);
        valid(inds) = false;
        [~,I] = max(ydata(zcr(inds)));
        valid(inds(I)) = true;
    end
end

zcr = zcr(valid);
        

inds = zcr;
