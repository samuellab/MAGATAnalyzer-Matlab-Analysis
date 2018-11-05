function [x,meany,sumy, numx] = meanyvsxFast (xdata, ydata, xaxis)
%function [x,meany, sumy, numx] = meanyvsxFast (xdata, ydata, xaxis)
%
%(xdata, ydata) form pairs, e.g. speed vs. angle
%
%for each interval in xaxis (xdata >= xaxis(j), <= xaxis(j+1)), 
%meany is the mean of all ydata with xdata in
%that interval
%

%meany = zeros([size(ydata,1) (length(xaxis)-1)]);
%sumy = meany;
%standarddeviation = meany;
%x = meany;

xdata = double(xdata); ydata = double(ydata);
if (size(ydata, 2) == 1)
    ydata = ydata';
end
if (size(xdata, 2) == 1)
    xdata = xdata';
end
[numx,bin] = histc(xdata, xaxis);
% %tic
% for j = 1:(length(xaxis) - 1)
% %    inds = find(xdata >= xaxis(j) & xdata < xaxis(j+1));
%     inds = bin == j;
%     sumy(:,j) = sum(ydata(:,inds),2);
%     x(j) = mean(xdata(inds));
%     meany(:,j) = mean(ydata(:,inds),2);
%     standarddeviation(:,j) = std(ydata(:,inds),0,2);
% end
% toc
% tic
meany = zeros([size(ydata,1) (length(xaxis)-1)]);


for k = 1:size(ydata, 1)
    meany(k,:) = accumarray(bin(bin > 0 & bin < length(xaxis))',ydata(k,bin > 0 & bin < length(xaxis))',[size(meany,2) 1],@mean)';
   
end
% A = accumarray(SUBS,VAL,SZ,FUN,FILLVAL)
%x = accumarray(bin(bin > 0 & bin < length(xaxis))',xdata(bin > 0 & bin < length(xaxis))',[size(meany,2) 1],@mean, NaN)';
%less accurate determination of bin center:
x = 0.5 * (xaxis(1:(end-1)) + xaxis(2:end));



numx = repmat(numx(1:(end-1)), size(ydata,1), 1);
sumy = meany.*numx;
