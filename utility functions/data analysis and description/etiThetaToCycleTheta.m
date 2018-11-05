function [theta, w, tx] = etiThetaToCycleTheta (theta, w, tx, tshift, period, trange)
%function [theta, w] = etiThetaToCycleTheta (theta, w, tx, tshift, period, trange)
%
%theta, w are mean and covariance estimates vs elapsed time; combines these
%into a mean and covariance estimate vs. cyclic time
%size theta is Ntime x ndim, size w is ndim x ndim x Ntime
deltaT = median(diff(tx));
if (~existsAndDefault('trange', []))
    trange(1) = tx(find(mod(tx + tshift, period) <= deltaT, 1, 'first'));
    trange = trange(1):period:max(tx);
    trange = [trange(1) trange(end)];
end
ti = tx >= min(trange) & tx <= max(trange);
tx = tx(ti);
tt = theta(ti,:);
ww = w(:,:,ti);
tc = mod(tx + tshift, period);
ii = round(tc / deltaT);
[~,tx] = meanyvsx(ii, tc, 0:(max(ii)+0.5));
theta = zeros(length(tx), size(tt,2));
w = zeros(size(ww,1), size(ww,2), length(tx));

n=1;
for j = 1:length(tx)
    
    ti = (ii + 1 == j);
    
    
    u = tt(ti, :)';
    wi = ww(:,:,ti);
    
    ni(j) = length(u);
    
    if(isempty(u) || isempty(wi))
        inds(n) = j;
        n=n+1;
        continue
    end
    
    for k = 1:size(wi,3)
        wi(:,:,k) = inv(wi(:,:,k));
        u(:,k) = wi(:,:,k)*u(:,k);
    end
    w(:,:,j) = inv(sum(wi,3));
    theta(j,:) = (w(:,:,j)*sum(u,2))';
end


return
% Stupid fix for the case when some bins have nothing in them-- please
% improve as you see fit -- Ruben 5/2017

if(n>1)
    for i=1:length(inds)
        if(inds(i)==1)
            theta(inds(i),:) = 0.5* ( theta(2,:) + theta(end,:) );
            w(:,:,inds(i)) = 0.5* (w(:,:,2) + w(:,:,end) );
        elseif(inds(i)==length(tx))
            theta(inds(i),:) = 0.5 * (theta(1,:) + theta(end-1,:));
            w(:,:,inds(i)) = 0.5 * ( w(:,:,1) + w(:,:,end-1) );
        else
            theta(inds(i),:) = 0.5* ( theta(inds(i)+1,:) + theta(inds(i)-1,:) );
            w(:,:,inds(i)) = 0.5* ( w(:,:,inds(i)+1) + w(:,:,inds(i)-1) );
        end
    end
end