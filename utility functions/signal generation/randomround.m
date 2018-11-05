function xr = randomround(x)

fx = floor(x);
dx = x - fx;
xr = fx + double(dx > rand(size(dx)));