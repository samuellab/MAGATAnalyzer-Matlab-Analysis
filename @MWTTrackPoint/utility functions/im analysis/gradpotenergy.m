function gpot = gradpotenergy (ctr, xderivim, yderivim)
%function gradpotenergy (ctr, xderivim, yderivim)
%interpolates along contour to find x derivatives and y derivatives.  
%if any part of contour is out of the image, gpot directs that point back
%towards the contour center

fx = interp2(xderivim, ctr(1,:), ctr(2,:), '*linear', NaN);
fy = interp2(yderivim, ctr(1,:), ctr(2,:), '*linear', NaN);

%nnz(~isfinite(fx))
fx(~isfinite(fx)) = mean(ctr(1,:)) - fx(~isfinite(fx));
fy(~isfinite(fy)) = mean(ctr(2,:)) - fy(~isfinite(fy));
gpot = [fx;fy];