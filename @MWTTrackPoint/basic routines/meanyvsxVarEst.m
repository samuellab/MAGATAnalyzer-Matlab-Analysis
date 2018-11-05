function [x,meany,v,numx] = meanyvsxVarEst (xdata, ydata, yvar, xaxis)
%function [x,meany,v,numx] = meanyvsxVarEst (xdata, ydata, yvar, xaxis)
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
if (length(xdata) ~= length(ydata))
    error ('data must be same length');
end

if (size(ydata, 1) > 1)
    if (size(ydata, 2) > 1)
        error ('multi d not written yet - time constraints');
    else
        ydata = ydata';
    end
end
if (size(yvar, 1) > 1)
    yvar = yvar';
end

meany = zeros([size(ydata,1) (length(xaxis)-1)]);
v = meany;

vi = 1./yvar;

yvi = vi.*ydata;

v(1,:) = accumarray(bin(bin > 0 & bin < length(xaxis))', vi(1,bin > 0 & bin < length(xaxis))',[size(meany,2) 1], @sum);
v = 1./v;
meany = accumarray(bin(bin > 0 & bin < length(xaxis))',yvi(1,bin > 0 & bin < length(xaxis))',[size(meany,2) 1],@sum)';
meany = v.*meany;
% 
% for k = 1:size(ydata, 1)
%     v(k,:) = accumarray
%     meany(k,:) = accumarray(bin(bin > 0 & bin < length(xaxis))',ydata(k,bin > 0 & bin < length(xaxis))',[size(meany,2) 1],@mean)';
%     if (nargout > 2)
%         standarddeviation(k,:) = accumarray(bin(bin > 0 & bin < length(xaxis))',ydata(k,bin > 0 & bin < length(xaxis))',[size(meany,2) 1],@std)';
%     end
% end
% A = accumarray(SUBS,VAL,SZ,FUN,FILLVAL)
x = accumarray(bin(bin > 0 & bin < length(xaxis))',xdata(bin > 0 & bin < length(xaxis))',[size(meany,2) 1],@mean, NaN)';
numx = repmat(numx(1:(end-1)), size(ydata,1), 1);
% sumy = meany.*numx;
% 
% if (nargout > 2)
%     standarderror = standarddeviation./sqrt(numx);
% end