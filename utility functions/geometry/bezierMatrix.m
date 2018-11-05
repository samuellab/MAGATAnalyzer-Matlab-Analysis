function [MS,dMS] = bezierMatrix(n, s)

MS = zeros(n+1, length(s));
a = s;
b = 1-s;
for k = 0:n
    MS(k+1,:) = nchoosek(n,k)*a.^k.*b.^(n-k);
end
if (nargout > 1)
    dMS = zeros(n+1, length(s));
    a = s;
    b = 1-s;
    dMS(1,:) = (-1)^(n)*((n)*b.^(n-1));
    dMS(n+1,:) = n*a.^(n-1);
    for k = 1:(n-1)
        dMS(k+1,:) = (-1)^(n-k)*nchoosek(n,k)*(k*a.^(k-1).*b.^(n-k) +(n-k)*a.^k.*b.^(n-k-1));
    end
end
