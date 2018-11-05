function x = makecolumn(x)
%function x = makecolumn(x)

if (size(x,1) == 1)
    x = x';
end