function [x,meany,standarderror,standarddeviation,sumy, numy] = meanyvsx_slidingwindow (xdata, ydata, xcenter, bindim, windowType, polar, varargin)
%function [x,meany,standarderror,standarddeviation, sumy, numy] = meanyvsx (xdata, ydata, xaxis, xcenter, bindim, windowType, polar, varargin)
%
%(xdata, ydata) form pairs, e.g. speed vs. angle
%
%for each point in xcenter, we take the weighted average of ydata with
%weighting determined by windowType and bindim
%
%if windowType is 'gaussian' (default), the weighting w_j is exp(-4*log(2)*(x_j -
%xc)^2/(bindim^2)) [bindim = FWHM]
%
%if windowType is 'step', the weighting w_j is 1 if abs(x_j - xc) < bindim/2
%and 0 otherwise [bindim = bin width]
%
%if windowType is 'halfgaussian', the weighting w_j is exp(-4*log(2)*(x_j -
%xc)^2/(bindim^2)) * (sign(x_j - xc) == sign(bindim)) :  in other words, if
%bindim > 0 we only consider xj > xc, and if bindim is < 0 we only consider
%xj < xc
%
% if windowType is 'none', we abort evaluation and return 0s for all values
%
%if polar is true, we add or sutbract multiples of 2*pi to the xdata to lie
%within +/-pi of the given center
%meany is the mean of all ydata with xdata in
%that interval
%
%passing 'period', period will add or subtract multiples of period to the
%xdata to lie within +/- period/2 of the given center:  polar = true  is
%the same as passing 'period', 2*pi
%
%
%standard error is the standard deviation of the data in the bin divided by
%the square root of the number of elements
%
%'clipNaN', true/[false] -- if true [x,meany,...] are clipped to
% the range for which isfinite(x) is true (aka, where some x data exists).
% this new option implemented 7/31/2012 may change behavior of other
% scripts

existsAndDefault('polar', false);
if (polar)
    period = 2*pi;
else
    period = [];
end
clipNaN = false;
varargin = assignApplicable(varargin);
sumy = zeros([size(ydata,1) length(xcenter)]);
meany = sumy;
numy = zeros(1, length(xcenter));
x = numy;
standarddeviation = meany;
standarderror = meany;

switch(lower(windowType))
    case 'gaussian'
        weightFun = @(xd,xc) exp(-4*log(2)*(xd-xc).^2./bindim^2);
        maxval = abs(2*bindim); %about 4.5 sigma
        minval = -abs(2*bindim);
    case 'halfgaussian'
        weightFun = @(xd,xc) exp(-4*log(2)*(xd-xc).^2./bindim^2).*(sign(xd-xc) == sign(bindim));
        maxval = 2*abs(bindim)*(bindim > 0);
        minval = -2*abs(bindim)*(bindim < 0);
    case 'step'
        weightFun = @(xd,xc) heavisideStep(bindim/2 - abs(xd-xc)); % We replace heaviside with heavisideStep so that the code can be compiled (B 2012/01/25)
        maxval = abs(bindim/2);
        minval = -abs(bindim/2);
    case 'none'
        return;
    otherwise
        disp('windowType must be gaussian or step');
        return;
        
end
valid = isfinite(xdata) & all(isfinite(ydata),1);
xdata = xdata(valid);
ydata = ydata(:,valid);

if (size(ydata, 1) > 1)

    for j = 1:length(xcenter)
        if (~isempty(period))
            xdata = mod(xdata - xcenter(j) + period/2, period) - period/2 + xcenter(j);
        end
        goodinds = (xdata - xcenter(j)) > minval & (xdata - xcenter(j)) < maxval;
        %goodinds = true(size(xdata));
        w = weightFun(xdata(goodinds),xcenter(j));

        numy(j) = sum(w);
        x(j) = sum(w.*xdata(goodinds))/numy(j);
       % w = repmat(w, size(ydata, 1), 1);
       for k = 1:size(ydata, 1)
            sumy(k,j) = sum(w.*ydata(k,goodinds), 2);
            meany(k,j) = sumy(k,j) / numy(j);
            standarddeviation(k,j) = sqrt(sum(w.*(ydata(k,goodinds) - meany(k,j)).^2)/numy(j));
            standarderror(k,j) =  standarddeviation(k,j)/sqrt(numy(j));
       end
    end
else
    for j = 1:length(xcenter)
        if (~isempty(period))
            xdata = mod(xdata - xcenter(j) + period/2, period) - period/2 + xcenter(j);
        end

        goodinds = (xdata - xcenter(j)) > minval & (xdata - xcenter(j)) < maxval;
       % goodinds = true(size(xdata));
        w = weightFun(xdata(goodinds),xcenter(j));

        numy(j) = sum(w);
        x(j) = sum(w.*xdata(goodinds))/numy(j);
        
        sumy(1,j) = sum(w.*ydata(goodinds));
        meany(1,j) = sumy(j) / numy(j);
        standarddeviation(:,j) = sqrt(sum(w.*(ydata(goodinds) - meany(1,j)).^2)/numy(j));
        standarderror(:,j) =  standarddeviation(:,j)/sqrt(numy(j));
    end
end

if (clipNaN)
    inds = isfinite(x);
    x = x(inds);
    meany = meany(:,inds);
    standarderror = standarderror(:,inds);
    standarddeviation = standarddeviation(:,inds);
    sumy = sumy(:,inds);
    numy = numy(inds);
end
