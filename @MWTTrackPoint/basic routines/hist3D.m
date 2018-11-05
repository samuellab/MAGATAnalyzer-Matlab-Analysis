function v = hist3D (xpts,ypts,zpts,xaxis,yaxis,zaxis)
%function v = hist3D (xpts,ypts,zpts,xaxis,yaxis,zaxis)
%
%makes a 3D histogram of xpts,ypts,zpts with bins defined by xaxis yaxis
%zaxis

v = zeros([length(yaxis),length(xaxis),length(zaxis)]);
inds = find (xpts >= min(xaxis) & xpts <= max(xaxis) & ypts >= min(yaxis) & ypts <= max(yaxis) & zpts >= min(zaxis) & zpts <= max(zaxis));
for n = inds
    j = find(xaxis < xpts(n),1,'last');
    k = find(yaxis < ypts(n),1,'last');
    l = find(zaxis < zpts(n),1,'last');
    v(k,j,l) = v(k,j,l)+1;
end

