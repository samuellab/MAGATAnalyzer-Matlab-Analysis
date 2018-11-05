function [midline, c1, c2] = splitOutline(cpts, endpt1, endpt2)
%function [midline, c1, c2] = splitOutline(cpts, enpt1, endpt2)
%
%splits contour into two halves, running from endpt1 to endpt2
%midline is the mean of these halves
%enpts can either be points, or indices into cpts
%cpts is 2Xnpts

cpts = double(cpts);
endpt1 = double(endpt1);
endpt2 = double(endpt2);


if numel(endpt1) == 1
   ind1 = endpt1;
else
   [~,ind1] = min(sum((cpts-repmat(endpt1, 1, size(cpts,2))).^2));
 %  ind1 = dt.nearestNeighbor(endpt1(1), endpt1(2));
end
if numel(endpt2) == 1
   ind2 = endpt2;
else
   [~,ind2] = min(sum((cpts-repmat(endpt2, 1, size(cpts,2))).^2));
   %ind2 = dt.nearestNeighbor(endpt2(1), endpt2(2));
end

if (ind1 < ind2)
    c1 = cpts(:,ind1:ind2);
    c2 = cpts(:,[ind1:-1:1 end:-1:ind2]);
else
    c1 = cpts(:, ind1:-1:ind2);
    c2 = cpts(:,[ind1:end 1:ind2]);
end

npts = max(length(c1), length(c2));
c1 = interp1(c1', linspace(1,length(c1), npts))';
c2 = interp1(c2', linspace(1,length(c2), npts))';
c1 = resampleContour(c1, 'closed', false);
c2 = resampleContour(c2, 'closed', false);
midline = resampleContour(0.5 * (c1 + c2), 'closed', false);
