function pctls = percentile (list, percentiles)
%function pctls = percentile (list, percentiles)
%
%finds the values corresponding to the percentiles (expressed as a fraction
%between (0 and 1) of list
%if list is 2D, list is flattened to 1D first

if (any (percentiles < 0 | percentiles > 1))
    disp ('percentiles must be expressed as a fraction between 0 and 1');
    pctls = [];
    return;
end
list = sort(list(:));
inds = round(percentiles*length(list));
inds(inds <= 0) = 1;
inds(inds > length(list)) = length(list);
pctls = list(inds);
