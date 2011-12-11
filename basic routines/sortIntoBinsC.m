function inds = sortIntoBinsC(Y, binEdges)
%function inds = sortIntoBinsC(Y, binEdges)
%
%inds(j) = k implies binEdges(k) <= Y(j) < binEdges(k+1)
%inds(j) = length(binEdges) implies Y(j) = binEdges(end)
%inds(j) = NaN implies Y(j) < binEdges(1) or Y(j) > binEdges(end)

inds = zeros(size(Y));
for j = 1:length(binEdges)
    inds = inds + (Y >= binEdges(j));
end

inds(inds == 0) = NaN;
inds(Y > binEdges(end)) = NaN;

