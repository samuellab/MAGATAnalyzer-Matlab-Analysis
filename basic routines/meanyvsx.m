function [x,meany,standarderror,standarddeviation,sumy] = meanyvsx (xdata, ydata, xaxis)
%function [x,meany,standarderror,standarddeviation, sumy] = meanyvsx (xdata, ydata, xaxis)
%
%(xdata, ydata) form pairs, e.g. speed vs. angle
%
%for each interval in xaxis (xdata >= xaxis(j), <= xaxis(j+1)), 
%meany is the mean of all ydata with xdata in
%that interval
%
%standard error is the standard deviation of the data in the bin divided by
%the square root of the number of elements
meany = zeros([size(ydata,1) (length(xaxis)-1)]);
sumy = meany;
for j = 1:(length(xaxis) - 1)
    inds = find(xdata >= xaxis(j) & xdata < xaxis(j+1));
    sumy(:,j) = sum(ydata(:,inds),2);
    x(j) = mean(xdata(inds));
    meany(:,j) = mean(ydata(:,inds),2);
    standarderror(:,j) = std(ydata(:,inds), 0, 2)/sqrt(length(inds));
    standarddeviation(:,j) = std(ydata(:,inds),0,2);
end