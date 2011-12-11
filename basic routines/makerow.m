function x = makerow(x)
%function x = makerow(x)

if (size(x,1) > 1)
    x = x';
end