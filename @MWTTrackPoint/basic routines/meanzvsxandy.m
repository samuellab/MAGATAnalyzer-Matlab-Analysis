function [meanz, numpts, sumz] = meanzvsxandy (xpts,ypts,zpts,xcenters,ycenters)
%function [meanz, numpts, sumz] = meanzvsxandy (xpts,ypts,zpts,xcenters,ycenters)

xcf = [xcenters 2*xcenters(end)-xcenters(end-1)];
xcr = [2*xcenters(1)-xcenters(2) xcenters];
xe = 0.5*(xcf + xcr);

ycf = [ycenters 2*ycenters(end)-ycenters(end-1)];
ycr = [2*ycenters(1)-ycenters(2) ycenters];
ye = 0.5*(ycf + ycr);

[~,jv] = histc(xpts, xe);
[~,iv] = histc(ypts, ye);

inds = (jv > 0 & jv < length(xe) & iv > 0 & iv < length(ye));

xpts = xpts(inds);
ypts = ypts(inds);
jv = jv(inds);
iv = iv(inds);
zpts = zpts(inds);

sumz = zeros(length(ycenters), length(xcenters));
numpts = sumz;
for j = 1:length(zpts)
    sumz(iv(j),jv(j)) = sumz(iv(j),jv(j)) + zpts(j);
    numpts(iv(j),jv(j)) = numpts(iv(j),jv(j)) + 1;
end
meanz = sumz./numpts;
meanz(numpts == 0) = 0;