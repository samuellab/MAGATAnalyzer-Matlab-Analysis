function inds = mylocalmax (vec)
%function inds = mylocalmax (vec)
%

inds = find(diff(vec) < 0 & ([0 diff(vec(1:(end-1)))] > 0));
