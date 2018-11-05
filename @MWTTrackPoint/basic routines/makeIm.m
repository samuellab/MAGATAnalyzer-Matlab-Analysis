function [im, xc, yc] = makeIm (xpts,ypts,xedges,yedges)
%function im = makeIm (xpts,ypts,xedges,yedges)
%
%makes a 2D histogram of xpts,ypts with bins defined by xedges yedges

%this code section modified
% from ndhist.m by  by Jonathan C. Lansey
% 30 Jan 2014 (Updated 15 Jul 2014) 
% 2D histogram which is actually kind of fast making use of matlab's histc

im = zeros(length(yedges),length(xedges));
[~,binX] = histc(xpts,xedges);
for ii=1:length(xedges)
    I = (binX==ii);
    N = histc(ypts(I),yedges);
    im(:,ii) = N';
end
%end snippet

im = im(1:(end-1), 1:(end-1));
xc = 0.5*(xedges(1:(end-1)) + xedges(2:end));
yc = 0.5*(yedges(1:(end-1)) + yedges(2:end));


%{
im = zeros([length(yaxis),length(xaxis)]);
inds = find (xpts >= min(xaxis) & xpts <= max(xaxis) & ypts >= min(yaxis) & ypts <= max(yaxis));
for n = inds'
    j = find(xaxis <= xpts(n),1,'last');
    k = find(yaxis <= ypts(n),1,'last');
    
    im(k,j) = im(k,j)+1;
end

%}