function [c2, lags ] = xc2(A,B,C, maxLag)
%function [c2, lags ] = xc2(A,B,C, maxLag)
%  c2(m,n) = sum_j (C(j)*A(lags(m)+j)*B(lags(n)+j)

nel = 2*maxLag + 1;
c2 = zeros(nel, nel);
lags = (-2*maxLag):maxLag;
for j = 1:nel
    m = -maxLag + j;
    if (m < 0)
        d = A(1:(end+m)).*C((-m + 1):end);
    else
        d = A((m+1):end).*C(1:(end-m));
    end
    c2(:,j) = xcorr(B,d, maxLag, 'none');
end


end

