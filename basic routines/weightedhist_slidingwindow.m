function [npts, weightedmean] = weightedhist_slidingwindow (xdata, xcenter, bindim, windowType, varargin)
%function [npts,weightedmean]  = meanyvsx (xdata, ydata, xaxis, xcenter, bindim, windowType, polar, varargin)
%
% returns the sum of the weight (npts)
% and the weighted mean sum(weight.*xdata)/npts
%
%if windowType is 'gaussian' (default), the weighting w_j is exp(-4*log(2)*(x_j -
%xc)^2/(bindim^2)) [bindim = FWHM]
%if windowType is 'step', the weighting w_j is 1 if abs(x_j - xc) < bindim/2
%and 0 otherwise [bindim = bin width]
%
%if windowType is 'halfgaussian', the weighting w_j is exp(-4*log(2)*(x_j -
%xc)^2/(bindim^2)) * (sign(x_j - xc) == sign(bindim)) :  in other words, if
%bindim > 0 we only consider xj > xc, and if bindim is < 0 we only consider
%xj < xc
%
%windowType may also be a function handle to the weighting function
%weight = windowType(xdata, xcenter, bindim)
%
%passing 'period', period will add or subtract multiples of period to the
%xdata to lie within +/- period/2 of the given center:  for polar data, pass 'period', 2*pi
%


period = [];

varargin = assignApplicable(varargin); %#ok<NASGU>

npts = zeros(1, length(xcenter));
weightedmean = npts;

if ischar(windowType)

    switch(lower(windowType))
        case 'gaussian'
            weightFun = @(xd,xc) exp(-4*log(2)*(xd-xc).^2./bindim^2);
        case 'halfgaussian'
            weightFun = @(xd,xc) exp(-4*log(2)*(xd-xc).^2./bindim^2).*(sign(xd-xc) == sign(bindim));
        case 'step'
            weightFun = @(xd,xc) heaviside(bindim/2 - abs(xd-xc));
        otherwise
            disp('windowType must be gaussian or step');
            return;

    end
else
    weightFun = @(xd, xc) windowType(xd, xc, bindim);
end

valid = isfinite(xdata);% & all(isfinite(ydata),1);
xdata = xdata(valid);

for j = 1:length(xcenter)
    if (~isempty(period))
        xdata = mod(xdata - xcenter(j) + period/2, period) - period/2 + xcenter(j);
    end
    
    w = weightFun(xdata,xcenter(j));
    npts(j) = sum(w);
    weightedmean(j) = sum(w.*xdata)/npts(j); 
    
end


