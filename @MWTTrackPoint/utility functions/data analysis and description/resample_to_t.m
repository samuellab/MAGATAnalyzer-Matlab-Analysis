function [y,ty] = resample_to_t (tx, x, tr)
%function [y,ty] = resample_to_t (tx, x, tr)
%uses the matlab resample command to map x onto tr
%takes care of end points; tr must be evenly spaced
tx = tx(:).';
tr = tr(:).';
trout = false;
if (size(x,2) ~= size(tx, 2))
    x = x.';
    trout = true;
end
if (size(x,2) ~= size(tx, 2))
    error ('length of x does not match length of tx');
end
xe = interp1(tx, x, tr([1 end]), 'nearest', 'extrap');
ii = tx > tr(1) & tx < tr(end);
txx = [tr(1) tx(ii) tr(end)];
xx = [xe(:,1) x(:,ii) xe(:,end)];
[y,ty] = resample(xx, txx, 1/mean(diff(tr)));
if (trout)
    y = y.';
end


end

