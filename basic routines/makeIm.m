function im = makeIm (xpts,ypts,xaxis,yaxis)
%function im = makeIm (xpts,ypts,xaxis,yaxis)
%
%makes a 2D histogram of xpts,ypts with bins defined by xaxis yaxis

im = zeros([length(yaxis),length(xaxis)]);
inds = find (xpts >= min(xaxis) & xpts <= max(xaxis) & ypts >= min(yaxis) & ypts <= max(yaxis));
for n = inds
    j = find(xaxis < xpts(n),1,'last');
    k = find(yaxis < ypts(n),1,'last');
    
    im(k,j) = im(k,j)+1;
end

