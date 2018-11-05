function centers = binCentersFromEdges (edges)
%function centers = binCentersFromEdges (edges)
%
%creates bin centers s.t. centers(j) is midway between (edges(j) and
%edges(j+1))
%

centers = 0.5 * (edges(1:end-1)+edges(2:end));