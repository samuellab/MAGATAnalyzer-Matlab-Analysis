function [theta_s,w_s ] = applySmoothingGaussianFilters( theta, w, q )
%function [theta_s,w_s ] = applySmoothingGaussianFilters( theta, w, q )
%
%applies smoothing step from Koyama, S., Eden, U.T., Brown, E.N., and Kass, R.E. (2009). Bayesian decoding of neural spike trains. Ann Inst Stat Math 62, 37.
%assuming F_k = identity matrix
%q is the diffusive noise term

theta_s = theta;
w_s = w;
for k = (length(theta_s)-1):-1:1
    if (size(q,3) > 1)
        qq = q(:,:,k);       
    else
        qq = q;
    end
    h = w(:,:,k)/(w(:,:,k) + qq);
    theta_s(:,k) = theta(:,k) + h * (theta_s(:,k+1) - theta(:,k));
    w_s(:,:,k) = w(:,:,k) + h * (w_s(:,:,k+1)-w(:,:,k)-qq)* h';
end

end

