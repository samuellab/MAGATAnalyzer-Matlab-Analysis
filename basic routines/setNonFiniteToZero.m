function x = setNonFiniteToZero(x)
%function x = setNonFiniteToZero(x)
%x(~isfinite(x)) = 0;

x(~isfinite(x)) = 0;