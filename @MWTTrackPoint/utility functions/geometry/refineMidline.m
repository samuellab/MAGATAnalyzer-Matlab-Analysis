function [midline, c1, c2] = refineMidline(midline, c1, c2)
%function refineMidline(midline, c1, c2)
%from an initial midline guess, finds an improved midline that 
%makes the midline perpendicular to the line between the contour halves

npts = length(midline);
midlinetangent = deriv(midline, npts/10);
ind1 = zeros([1 npts]);
ind2 = ind1;
for j = 1:npts
    d1 = c1 - repmat(midline(:,j), [1 npts]);
    d2 = c2 - repmat(midline(:,j), [1 npts]);
    mt = repmat(midlinetangent(:,j), [1 npts]);
    [~,I] = min(abs(sum(d1.*mt)));
    ind1(j) = I;
    [~,I] = min(abs(sum(d2.*mt)));
    ind2(j) = I;
end
c1 = c1(:,ind1);
c2 = c2(:,ind2);
midline = 0.5 * (c1 + c2);