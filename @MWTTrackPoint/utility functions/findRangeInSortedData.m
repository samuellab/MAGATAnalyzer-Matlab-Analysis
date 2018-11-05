function inds = findRangeInSortedData(x, xmin, xmax)
%find all values of x in range [xmin, xmax)

a = 1;
b = length(x);
if (x(end) < xmin)
    inds = [];
    return;
end
if (x(1) >= xmax)
    inds = [];
    return;
end

while (b - a > 1)
    ii = floor ((a+b)/2);
    if (x(ii) < xmin)
        a = ii;
    else
        b = ii;
    end
end
if (x(a) >= xmin) %handle edge case
    ind1 = a;
else
    ind1 = b;
end
if (x(ind1) >= xmax)
    inds = [];
    return;
end
a = ind1;
b = length(x);
while (b - a > 1)
    ii = ceil ((a+b)/2);
    if (x(ii) >= xmax)
        b = ii;
    else
        a = ii;
    end
end
if (x(b) < xmax) %handle edge case
    ind2 = b;
else
    ind2 = a;
end
inds = ind1:ind2;

