function edges = binEdgesFromCenters (centers)
%function edges = binEdgesFromCenters (centers)
%
%creates bin edges s.t. centers(j) is midway between (edges(j) and
%edges(j+1))
%
%modified from MATLAB's built in hist code

%matlabs code, which doesn't actually produce bins centered between edges,
%if bin centers are non uniform
 xx = centers(:)';
 binwidth = diff(xx); binwidth(end + 1) = binwidth(end);
 xx = [xx(1)-binwidth(1)/2 xx+binwidth/2];
% xx(1) = -Inf;
% xx(end) = Inf;
 edges = xx;