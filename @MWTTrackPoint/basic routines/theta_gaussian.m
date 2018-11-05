function g = theta_gaussian (x, m, s)
%function g = theta_gaussian (x, m, s)
%calculates gaussian using theta math

td = mod(x - m, 2*pi);

td = min(abs(td), abs(2*pi - td));

g = (1/sqrt(2*pi*s))*exp(-td.^2/(2*s));
dx = diff(x);
dx(end+1) = dx(end);
g = g.*dx;
