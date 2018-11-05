function [x,meany,standarderror,standarddeviation,sumy, numx] = meanyvsx (xdata, ydata, xaxis)
%function [x,meany,standarderror,standarddeviation, sumy, numx] = meanyvsx (xdata, ydata, xaxis)
%
%(xdata, ydata) form pairs, e.g. speed vs. angle
%
%for each interval in xaxis (xdata >= xaxis(j), <= xaxis(j+1)), 
%meany is the mean of all ydata with xdata in
%that interval
%
%standard error is the standard deviation of the data in the bin divided by
%the square root of the number of elements

%meany = zeros([size(ydata,1) (length(xaxis)-1)]);
%sumy = meany;
%standarddeviation = meany;
%x = meany;

xdata = double(xdata); ydata = double(ydata);

if (size(xdata, 2) == 1)
    xdata = xdata';
end
if (size(ydata, 2) ~= size(xdata,2))
    ydata = ydata';
end
[numx,bin] = histc(xdata, xaxis);

if (length(xdata) ~= length(ydata))
    error ('data must be same length');
end

meany = zeros([size(ydata,1) (length(xaxis)-1)]);

standarddeviation = meany;

for k = 1:size(ydata, 1)
    meany(k,:) = accumarray(bin(bin > 0 & bin < length(xaxis))',ydata(k,bin > 0 & bin < length(xaxis))',[size(meany,2) 1],@mean)';
    if (nargout > 2)
        standarddeviation(k,:) = accumarray(bin(bin > 0 & bin < length(xaxis))',ydata(k,bin > 0 & bin < length(xaxis))',[size(meany,2) 1],@std)';
    end
end
% A = accumarray(SUBS,VAL,SZ,FUN,FILLVAL)
x = accumarray(bin(bin > 0 & bin < length(xaxis))',xdata(bin > 0 & bin < length(xaxis))',[size(meany,2) 1],@mean, NaN)';
numx = repmat(numx(1:(end-1)), size(ydata,1), 1);
sumy = meany.*numx;

if (nargout > 2)
    standarderror = standarddeviation./sqrt(numx);
end