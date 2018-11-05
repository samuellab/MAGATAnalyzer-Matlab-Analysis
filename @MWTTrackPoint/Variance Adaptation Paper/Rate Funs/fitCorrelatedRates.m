function [ fP, LL, rates, fitRates, varHL] = fitCorrelatedRates( varstruct, opstruct, degree, coordinates )
%function [ fP, LL, rates, fitRates, varHL] = fitCorrelatedRates( varstruct, opstruct, degree, coordinates )

% % outputs:
% fP are fit params for different rate models
% LL are log-likelihoods for different rate models
% rates are the actual rate functions
% fitRates are the fit functions for different models
% varHL are the var structures containing turn and run data

% inputs:
% degree = 1 or 2 for rate functions being exponentials of a 1st or 2nd degree polynomial
% coordinates: 'OL' fits different rate functions in odor-light space
%              'UV' fits different rate functions in u-v space (u=cos(th)*xo+sin(th)*xl)



theta0 = deg2rad(45); % initial guess

op = optimoptions('fminunc');
op.Display = 'off';
op.Algorithm = 'quasi-newton';

nBins = opstruct.numLxBins;
deltaT = median(diff(varstruct(1).fullensemble.eti));


for i=1:length(varstruct)
    btdvar = varstruct(i);
    tdata(i,:) = btdvar.turn.x_conv;
    rdata(i,:) = btdvar.noturn.x_conv;
    tton = btdvar.turn.ton;
    ttoff = btdvar.turn.toff;
    teti = btdvar.turn.eti;
    rton = btdvar.noturn.ton;
    rtoff = btdvar.noturn.toff;
    reti = btdvar.noturn.eti;
end

tH = tton > opstruct.adaptationTime & tton < ttoff & teti >= min(opstruct.timeRange) & teti < max(opstruct.timeRange);
tL = ttoff > opstruct.adaptationTime & ttoff < tton & teti >= min(opstruct.timeRange) & teti < max(opstruct.timeRange);
tU = tton > 0 & tton < opstruct.adaptationTime & teti >= min(opstruct.timeRange) & teti < max(opstruct.timeRange);
tD = ttoff > 0  & ttoff < opstruct.adaptationTime & teti >= min(opstruct.timeRange) & teti < max(opstruct.timeRange);
rH = rton > opstruct.adaptationTime & rton < rtoff & reti >= min(opstruct.timeRange) & reti < max(opstruct.timeRange);
rL = rtoff > opstruct.adaptationTime & rtoff < rton & reti >= min(opstruct.timeRange) & reti < max(opstruct.timeRange);
rU = rton > 0 & rton < opstruct.adaptationTime & reti >= min(opstruct.timeRange) & reti < max(opstruct.timeRange);
rD = rtoff > 0  & rtoff < opstruct.adaptationTime & reti >= min(opstruct.timeRange) & reti < max(opstruct.timeRange);

for i=1:length(varstruct)
    tdatH(i,:) = tdata(i,tH);
    tdatL(i,:) = tdata(i,tL);
    rdatH(i,:) = rdata(i,rH);
    rdatL(i,:) = rdata(i,rL);
end

if(degree==1)
    p0 = @(xt, xr) [mean(xt)/var(xr) log(length(xt)/(length(xr)*deltaT))];
else
    p0 = @(xt, xr) [-0.05 mean(xt)/var(xr) log(length(xt)/(length(xr)*deltaT))];
end

if(strcmpi(coordinates, 'OL'))
    
    rate = @(rateP, dat) exp(polyval(rateP, dat));
    rate_resc = @(a, rateP, dat) exp(polyval(rateP, a.*dat));
    rate2d = @(rP, dat) exp( rP(1) * dat(1,:) + rP(2) * dat(2,:) + rP(3) );
    
    nlogP = @(p, xt, xr) -sum(log(rate(p, xt) )) + sum(rate(p, xr)*deltaT);
    nlogP2d = @(p, xt, xr) -sum(log(rate2d(p, xt))) + sum(rate2d(p, xr)*deltaT);
    
    poH = p0(tdatH(1,:), rdatH(1,:));
    plH = p0(tdatH(2,:), rdatH(2,:));
    poL = p0(tdatL(1,:), rdatL(1,:));
    plL = p0(tdatL(2,:), rdatL(2,:));
    
    % rH(o,l) = boH * xoH + blH * xlH + c;
    ff = @(p) nlogP2d(p, tdatH, rdatH); % p = [boH, blH, c]
    ff0 = [poH(1) plH];
    tic
    [fP.olH, LL.olH] = fminunc(ff, ff0, op);
    toc
    
    % rL(o,l) = boL * xoL + blL * xlH + c;
    ff = @(p) nlogP2d(p, tdatL, rdatL); %  p = [boL, blL, c]
    ff0 = [poL(1) plL];
    ff0 = fP.olH;
    tic
    [fP.olL, LL.olL] = fminunc(ff, ff0, op);
    toc
    
    % rH(o,l) = boH * xoH + blH * xlH + c; rL(o,l) = boL * xoL + blL * xlH + c;
    ff = @(p) nlogP2d( [p(1:2) p(end)], tdatH, rdatH) + nlogP2d( p(3:end), tdatL, rdatL); % p = [boH, blH, boL, blL, c]
    ff0 = [fP.olH(1:2) fP.olL];
    tic
    [fP.ol, LL.ol] = fminunc(ff, ff0, op);
    toc
    
    
    
elseif(strcmpi(coordinates, 'UV'))
    UVdat = @(xdat, theta) [cos(theta)*xdat(1,:) + sin(theta)*xdat(2,:); -sin(theta)*xdat(1,:) + cos(theta)*xdat(2,:)];
    Udat = @(xdat, theta) cos(theta)*xdat(1,:) + sin(theta)*xdat(2,:);
    Vdat = @(xdat, theta) -sin(theta)*xdat(1,:) + cos(theta)*xdat(2,:);
    
    rate = @(rateP, dat) exp(polyval(rateP, dat));
    rate_resc = @(a, rateP, dat) exp(polyval(rateP, a.*dat));
    ratefunUV = @(theta, rPu, rPv, tdat) exp(polyval(rPu, Udat(tdat, theta))).*exp(rPv*Vdat(tdat, theta));
    ratefunUV_resc = @(aU, aV, theta, rPu, rPv, dat) exp(polyval(rPu, aU.*Udat(dat, theta))).*exp(rPv*aV.*Vdat(dat, theta));
    
    nlogPU = @(p, xt, xr) -sum(log(rate(p(2:end), Udat(xt, p(1)) ) )) + sum(rate(p(2:end), Udat(xr, p(1)))*deltaT);
    nlogPV = @(p, xt, xr) -sum(log(rate(p(2:end), Vdat(xt, p(1)) ) )) + sum(rate(p(2:end), Vdat(xr, p(1)))*deltaT);
    nlogPU_resc = @(p, xt, xr) -sum(log(rate_resc(p(1), p(3:end), Udat(xt, p(2)) ) )) + sum(rate_resc(p(1), p(3:end), Udat(xr, p(2)))*deltaT);
    nlogPV_resc = @(p, xt, xr) -sum(log(rate_resc(p(1), p(3:end), Vdat(xt, p(2)) ) )) + sum(rate_resc(p(1), p(3:end), Vdat(xr, p(2)))*deltaT);
    nlogPuv = @(p, xt, xr) -sum(log(ratefunUV(p(1), p(2:end-1), p(end), xt))) + sum(ratefunUV(p(1), p(2:end-1), p(end), xr)*deltaT);
    nlogPuv_resc = @(p, xt, xr) -sum(log(ratefunUV_resc(p(1), p(2), p(3), p(4:end-1), p(end), xt))) + sum(ratefunUV_resc(p(1), p(2), p(3), p(4:end-1), p(end), xr)*deltaT);
    
    
    tUH = Udat(tdatH, theta0);
    rUH = Udat(rdatH, theta0);
    tUL = Udat(tdatL, theta0);
    rUL = Udat(rdatL, theta0);
    tVH = Vdat(tdatH, theta0);
    rVH = Vdat(rdatH, theta0);
    
    
    puH = p0(tUH, rUH);
    pvH = p0(tVH, rVH);
    puL = p0(tUL, rUL);
    puAll = p0([tUH tUL], [rUH rUL]);
    
    % 1: r(uH) = a*u_H(theta) + b
    tic
    ff = @(p) nlogPU(p, tdatH, rdatH);  % p = [theta, a, b]
    ff0 = [theta0 puH];
    [fP.uH, LL.uH] = fminunc(ff, ff0, op);
    toc
    
    theta = fP.uH(1);
    rates.u_H = UVrates(tdatH, rdatH, theta, deltaT, nBins, 1);
    fitRates.u_H = rate( fP.uH(2:end), rates.u_H.lx );
    
    % 2: r(uL) = a*u_L(theta) + b
    tic
    ff = @(p) nlogPU(p, tdatL, rdatL);  % p = [theta, a, b]
    ff0 = [theta0 puL];
    [fP.uL, LL.uL] = fminunc(ff, ff0, op);
    toc
    theta = fP.uL(1);
    rates.u_L = UVrates(tdatL, rdatL, theta, deltaT, nBins, 1);
    fitRates.u_L = rate( fP.uL(2:end), rates.u_L.lx );
    
    % 3: fit theta_low with fixed rateP(high)
    tic
    ff = @(p) nlogPU([p, fP.uH(2:end)], tdatL, rdatL);  % p = [theta, a, b]
    ff0 = [theta0];
    [fP.thetaL, LL.uL_fixedrateH] = fminunc(ff, ff0, op);
    toc
    
    
    % 4: r(uH) = a * uH(theta) + b; r(uL) = a * uL(theta) + b;
    ff = @(p) nlogPU(p, tdatH, rdatH) + nlogPU(p, tdatL, rdatL); % p = [theta, a, b]
    ff0 = [fP.uH];
    tic
    [fP.u_noresc, LL.u_noresc] = fminunc(ff, ff0, op);
    toc

    % 5: r(uH) = a * uH(thetaH) + b; r(uL) = a * uL(thetaL) + b;
%     ff = @(p) nlogPU([p(1) p(3:end)], tdatH, rdatH) + nlogPU([p(2) p(3:end)], tdatL, rdatL);
%     ff0 = [theta0 theta0 fP.u_noresc(2:end)];
%     tic
%     [fP.u_theta, LL.u_theta] = fminunc(ff, ff0, op); % p = [thetaH, thetaL, a, b]
%     toc
    
    % r(uH) = aH * uH(theta) + b; r(uL) = aL * uL(theta)) + b; - same as resc. model below
    % ff = @(p) nlogPU([p(1) p(2) p(end)], tdatH, rdatH) + nlogPU([p(1) p(3) p(end)], tdatL, rdatL);
    % ff0 = [theta0 1 1 fP.u_noresc(end)];
    % tic
    % [fP.u_lin, LL.u_lin] = fminunc(ff, ff0, op); % p = [theta, aH, aL, b]
    % toc
    
    % 5: r(u) = lambda(alpha_u * u)
    % r(uH) = a * uH(theta) + b; r(uL) = a * ( alphaU * uL(theta)) + b;
    ff = @(p) nlogPU(p(2:end), tdatH, rdatH) + nlogPU_resc(p, tdatL, rdatL);
    ff0 = [1 fP.u_noresc];
    tic
    [fP.u_resc, LL.u_resc] = fminunc(ff, ff0, op); % params = [alphaU, theta, a, b]
    toc
    theta = fP.u_resc(2);
    rates.uH_resc = UVrates(tdatH, rdatH, theta, deltaT, nBins, 1);
    rates.uL_resc = UVrates(tdatL, rdatL, theta, deltaT, nBins, 1);
    fitRates.uH_resc = rate_resc( 1, fP.u_resc(3:end), rates.uH_resc.lx );
    fitRates.uL_resc = rate_resc( fP.u_resc(1), fP.u_resc(3:end), rates.uL_resc.lx );
    
    % 6: r(vH) = a * vH(theta) + b; r(vL) = a * ( alphaV * vL(theta)) + b;
    ff = @(p) nlogPV(p(2:end), tdatH, rdatH) + nlogPV_resc(p, tdatL, rdatL);
    ff0 = [1 theta0 pvH];
    tic
    [fP.v_resc, LL.v_resc] = fminunc(ff, ff0, op); % p = [alphaV, theta, a, b]
    toc
    
    % r(uH) = a * uH(thetaH) + b; r(uL) = a * ( alpha_L * uL(thetaL)) + b;
%     ff = @(p) nlogPU([p(1) p(4:end)], tdatH, rdatH) + nlogPU_resc(p(2:end), tdatL, rdatL);
%     ff0 = [theta0 fP.u_resc];
%     tic
%     [fP.u_theta_resc, LL.u_theta_resc] = fminunc(ff, ff0, op); % p = [thetaH, alphaU, thetaL, a, b]
%     toc
 

% 7: r(u,v) = lambda(u, v) with fixed theta:
% r(uH,vH) = a * u_H(theta) + b + c * v_H(theta);  r(uL,vL) = a * u_L(theta) + b + c * v_L(theta)
theta = fP.u_resc(2);
ff = @(p) nlogPuv_resc([1 1 theta p(1:end)], tdatH, rdatH) + nlogPuv_resc([1 1 theta p(1:end)], tdatL, rdatL);
ff0 = [fP.u_resc(3:end) fP.v_resc(3:end-1)];
tic
[fP.uv, LL.uv] = fminunc(ff, ff0, op); % params = [a, b, c] 
toc

% 8: r(u,v) = lambda(alpha_u*u, alpha_v*v) with fixed theta:
% r(uH,vH) = a * u_H(theta) + b * v_H(theta) + c ;  r(uL,vL) = alpha_u *(a * u_L(theta)) + alpha_v *(b * v_L(theta))+ c
theta = fP.u_resc(2);
ff = @(p) nlogPU([theta p(3:3+degree-1) p(end)], tdatH, rdatH) + nlogPU_resc([p(1) theta p(3:3+degree-1) p(end)], tdatL, rdatL) + nlogPV([theta p(3+degree:3+2*degree-1) p(end)], tdatH, rdatH) + nlogPV_resc([p(2) theta p(3+degree:3+2*degree-1) p(end)], tdatL, rdatL);
ff0 = [fP.u_resc(1) 1 fP.u_resc(3:end-1) fP.v_resc(3:end)];
tic
[fP.uv_resc_fixedTheta, LL.uv_resc_fixedTheta] = fminunc(ff, ff0, op); % params = [alpha_u, alpha_v, a, b, c] 
toc

end
end

function rf = UVrates(tdat, rdat, theta, deltaT, nBins, stimNum)

if(~isempty(stimNum))
    if(stimNum == 1)
        NewDat = @(xdat, theta) cos(theta)*xdat(1,:) + sin(theta)*xdat(2,:);
    elseif(stimNum==2)
        NewDat = @(xdat, theta) -sin(theta)*xdat(1,:) + cos(theta)*xdat(2,:);
    end
    
    xt = NewDat(tdat, theta);
    xr = NewDat(rdat, theta);
    lx = linspace(-3*std(xr), 3*std(xr), nBins);
    
    
    nt = histcounts(xt, binEdgesFromCenters(lx));
    nr = histcounts(xr, binEdgesFromCenters(lx));
    
    rf.lx = lx;
    rf.rate = 1/deltaT * nt ./ nr;
    rf.rate_eb = 1/deltaT * sqrt(nt) ./ nr;
    
else
    UVdat = @(xdat, theta) [cos(theta)*xdat(1,:) + sin(theta)*xdat(2,:); -sin(theta)*xdat(1,:) + cos(theta)*xdat(2,:)];
    
    xt = UVdat(tdat, theta);
    xr = UVdat(rdat, theta);
    td(:,1) = xt(1,:);
    td(:,2) = xt(2,:);
    rd(:,1) = xr(1,:);
    rd(:,2) = xr(2,:);
    lx(:,1) = linspace(-3*std(rd(:,1)), 3*std(rd(:,1)), nBins);
    lx(:,2) = linspace(-3*std(rd(:,2)), 3*std(rd(:,2)), nBins);
    [lxx, lxy] = meshgrid(lx(:,1), lx(:,2));
    lxdata = [lxx(:) lxy(:)];
    
    
    
    h = makeIm(rd(:,1), rd(:,2), binEdgesFromCenters(lx(:,1)), binEdgesFromCenters(lx(:,2)));
    ht = makeIm(td(:,1), td(:,2), binEdgesFromCenters(lx(:,1)), binEdgesFromCenters(lx(:,2)));
    rhist = reshape(h, size(lxx));
    thist = reshape(ht, size(lxx));
    
    rf.lx = lx';
    rf.rate(1,:) = 1/deltaT * sum(thist, 1)./sum(rhist, 1);
    rf.rate_eb(1,:) = 1/deltaT * sqrt(sum(thist, 1))./sum(rhist, 1);
    rf.rate(2,:) = 1/deltaT * sum(thist, 2)./sum(rhist, 2);
    rf.rate_eb(2,:) = 1/deltaT * sqrt(sum(thist, 2))./sum(rhist, 2);
end
end
