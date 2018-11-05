function [x,mediany] = medianyvsx (xdata, ydata, xaxis)
%function [x,mediany] = medianyvsx (xdata, ydata, xaxis)
%
%(xdata, ydata) form pairs, e.g. speed vs. angle
%
%for each interval in xaxis, mediany is the median of all ydata with xdata in
%that interval

for j = 1:(length(xaxis) - 1)
    inds = find(xdata >= xaxis(j) & xdata < xaxis(j+1));
    x(j) = mean(xdata(inds));
    mediany(j) = median(ydata(inds));
end