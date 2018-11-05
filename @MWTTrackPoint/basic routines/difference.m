function dv = difference (vector, spacing)
%function dv = difference (vector, spacing)
%
%dv(j) = vector(j+spacing)-vector(j) if j <= length(vector) - spacing
%otherwise dv(j) = (vector(end)-vector(j))*spacing/(length(vector) - j)
%dv(end) == dv(end-1);
dv = zeros(size(vector));
n = length(vector) - spacing;
dv(1:n) = vector(1+spacing:end) - vector(1:n);
for k = (n+1):(length(vector)-1)
    dv(k) = (vector(end) - vector(k))*spacing/(length(vector)-k);
end
dv(end) = dv(end-1);
    