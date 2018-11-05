function saw = mysawtooth(x, a, low, high)
%function saw(x)
%
%produces a saw wave with angular frequency 1 and upward duty cycle a
% optional args, low and high provide the min and max values
if (nargin < 2)
    a = 1;
end
if (nargin < 3)
    low = -1;
end
if (nargin < 4)
    high = -low;
end

saw = mod(x, 2*pi)/(2*pi);
up = logical(saw <= a+eps);
saw(up) = saw(up)/a;
saw(~up) = (1 -saw(~up))/(1 - a);

saw = saw*(high - low) + low;

