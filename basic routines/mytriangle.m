function tri = mytriangle(x)
%function mytriangle(x)
%
%produces a triangle wave with amplitude 1 and angular frequency 1
xi = linspace(x(1)-2*pi, x(end)+2*pi, 10*length(x)+100);
square = sign(sin(xi));
inds = find(diff(square) ~= 0);
tri = interp1(xi(inds), square(inds), x, 'linear');