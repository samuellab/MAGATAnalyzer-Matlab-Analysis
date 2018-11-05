function [wholeinds, fracinds] = findPositiveZeroCrossings(x)
%function [wholeinds, fracinds] = findPositiveZeroCrossings(x)
%
%finds locations where x(i) < 0 and x(i+1) > 0;  
%wholeinds -- i
%fracinds -- fractional index where x(fracind) = 0 by linear interpolation
%if x(i) == 0, then there is a positive zero crossing iff
%  x(i - 1) < 0 and
%  x(i + 1) > 0

s = sign(x);
if (all (s == 0))
    wholeinds = [];
    fracinds = [];
    return;
end
if (s(end) == 0)
    s(end) = s(find(s ~= 0, 1, 'last'));
end
while any(s == 0)
   s(s == 0) = s(find(s == 0) + 1);
end

wholeinds = find(diff(s) > 0);
a = x(wholeinds);
b = x(wholeinds + 1);
fracinds = wholeinds - a./(b-a);
