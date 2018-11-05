function [rate, xaxis, optw, rate_ci] = rateAdaptiveBW (eventdata, normdata, xaxis, kernelType, causal, alpha)
% function [rate, rate_ci] = rateAdaptiveBW (eventdata, normdata, xaxis, kernelType, causal, alpha)
% eventdata = xvalue of events (e.g. time of turns, heading angle at time
% of turn)
% normdata = xvalue when event could have occurred (e.g. time animal was in
% a run)
% xaxis = axis for making histograms (bin centers)
% rate is in terms of the time bin size at which normdata was acquired
% e.g. a rate of 1 means 1 event per time bin - usually you will want to
% multiply this rate by the sampling rate to get the rate in hertz

existsAndDefault('xaxis', []);
existsAndDefault('kernelType', 'Gauss');
existsAndDefault('causal', false);
existsAndDefault('alpha', 0.95);


if (isempty(xaxis))
    [~,xaxis, optw] = ssvkernel(eventdata);
else 
    xaxis = linspace(min(xaxis), max(xaxis), length(xaxis));
    [~,t,optw] = ssvkernel(eventdata, xaxis);
    optw = interp1(t, optw, xaxis, 'linear', 'extrap');
end
if (nargout > 3)
    [h1, h1_ci] =  histVarBW (eventdata, xaxis, optw, kernelType, causal, alpha);
else
    h1 =  histVarBW (eventdata, xaxis, optw, kernelType, causal, alpha);
end
h2 =  histVarBW (normdata, xaxis, optw, kernelType, causal, alpha);
rate = h1 ./ h2;
if (nargout > 3)
    rate_ci(1,:) = h1_ci(1,:)./h2;
    rate_ci(2,:) = h1_ci(2,:)./h2;
end
