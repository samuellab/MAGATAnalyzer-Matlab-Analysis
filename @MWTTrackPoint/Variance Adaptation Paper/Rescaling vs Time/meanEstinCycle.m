function [alpha_ton, valpha_ton, tx_ton] = meanEstinCycle (alpha, valpha, tx, tshift, period, trange, timeType)
% function [alpha_ton, valpha_ton, tx_ton] = meanEstinCycle (alpha, valpha, tx, tshift, period, trange, timeType)

% if timeType='eti' : calculates mean + error vs. time in cycle (t=0:period) from estimates+errors vs. eti (e.g. t=60:1200)
% if timeType='ton' : averages mean + error vs. time in cycle (t=0:period)



deltaT = median(diff(tx));

if(size(alpha,2)>size(alpha,1))
    alpha = alpha';
end
nstim = size(alpha,2);

if(strcmpi(timeType, 'ton'))
    trange = [min(tx) max(tx)];
    tx = repmat(tx, [1 round(size(alpha, 1)/length(tx))]);
end

if (~existsAndDefault('trange', []))
    trange(1) = tx(find(mod(tx + tshift, period) <= deltaT, 1, 'first'));
    trange = trange(1):period:max(tx);
    trange = [trange(1) trange(end)];
end
ti = tx >= min(trange) & tx <= max(trange);
tx = tx(ti);
tt = alpha(ti,:);

if(strcmpi(timeType, 'ton'))
    tc = tx;
else
    tc = mod(tx + tshift, period);
end
ii = round(tc / deltaT);
[~,tx_ton] = meanyvsx(ii, tc, 0:(max(ii)+0.5));
alpha_ton = zeros(length(tx_ton), size(tt,2));
% valpha_ton = zeros(size(valpha));

if(nstim==2)
    
    ww = valpha(:,:,ti);
    valpha_ton = zeros(size(ww,1), size(ww,2), length(tx_ton));
    
    n=1;
    for j = 1:length(tx_ton)
        
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
        valpha_ton(:,:,j) = inv(sum(wi,3));
        alpha_ton(j,:) = (valpha_ton(:,:,j)*sum(u,2))';
    end
    alpha_ton = alpha_ton';
    
elseif(nstim==1)
    
    valpha = squeeze(valpha);
    ww = valpha(ti);
    valpha_ton = zeros([1,length(tx_ton)]);
    
    n=1;
    for j = 1:length(tx_ton)
        
        ti = (ii + 1 == j);
        
        
        u = tt(ti, :)';
        wi = ww(ti);
        
        ni(j) = length(u);
        
        if(isempty(u) || isempty(wi))
            inds(n) = j;
            n=n+1;
            continue
        end
        
        for k = 1:length(wi)
            wi(k) = 1/wi(k);
            u(:,k) = wi(k)*u(:,k);
        end
        inds = ~isinf(wi);
        valpha_ton(j) = inv(sum(wi(inds)));
        alpha_ton(j,:) = (valpha_ton(j)*sum(u(inds),2))';
    end
    alpha_ton = alpha_ton';
    
    
    
end
