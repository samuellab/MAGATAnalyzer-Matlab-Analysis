function [ fP, LL ] = fitLowHighRates( lowvar, highvar, opstruct, degree, rescale_type )
%function [ fP, LL ] = fitLowHighRates( lowvar, highvar, opstruct, degree, rescale_type )

% fit one-stimulus data to rescaling rate function
% rescale_type = {'input', 'output', 'output_rect'} or any combination of those
% degree: rate function polynomial degree in the exponential: degree=1: r=exp(ax+b) or degree=2: r=exp(ax^2+bx+c)

op = optimoptions('fminunc');
op.Display = 'off';
op.Algorithm = 'quasi-newton';

nBins = opstruct.numLxBins;
deltaT = median(diff(lowvar.var.fullensemble.eti));

tetiH = highvar.var.turn.eti;
tetiL = lowvar.var.turn.eti;
retiH = highvar.var.noturn.eti;
retiL = lowvar.var.noturn.eti;

tH = tetiH >= min(opstruct.timeRange) & tetiH < max(opstruct.timeRange);
tL = tetiL >= min(opstruct.timeRange) & tetiL < max(opstruct.timeRange);
rH = retiH >= min(opstruct.timeRange) & retiH < max(opstruct.timeRange);
rL = retiL >= min(opstruct.timeRange) & retiL < max(opstruct.timeRange);

tdatH = highvar.var.turn.x_conv(tH);
tdatL = lowvar.var.turn.x_conv(tL);
rdatH = highvar.var.noturn.x_conv(rH);
rdatL = lowvar.var.noturn.x_conv(rL);

rates.lxL = linspace(-3*std(rdatL), 3*std(rdatL), nBins);
rates.lxH = linspace(-3*std(rdatH), 3*std(rdatH), nBins);
rates.lx_all = linspace(-3*std([rdatL; rdatH]), 3*std([rdatL; rdatH]), nBins);
ntL = histcounts(tdatL, binEdgesFromCenters(rates.lxL));
nrL = histcounts(rdatL, binEdgesFromCenters(rates.lxL));
ntH = histcounts(tdatH, binEdgesFromCenters(rates.lxH));
nrH = histcounts(rdatH, binEdgesFromCenters(rates.lxH));
nt_all = histcounts([tdatL; tdatH], binEdgesFromCenters(rates.lx_all));
nr_all = histcounts([rdatL; rdatH], binEdgesFromCenters(rates.lx_all));
rates.rL = 1./deltaT *ntL./nrL;
rates.rH = 1./deltaT *ntH./nrH;
rates.rall = 1./deltaT*nt_all./nr_all;


if(degree==1)
    p0 = @(xt, xr) [mean(xt)/var(xr) log(length(xt)/(length(xr)*deltaT))];
else
    p0 = @(xt, xr) [-0.05 mean(xt)/var(xr) log(length(xt)/(length(xr)*deltaT))];
end




    
    for k=1:length(rescale_type)
        
        if(strcmpi(rescale_type{k}, 'input'))
            ratefun = @(alpha, rateP, dat) exp(polyval(rateP, alpha.*dat));
        elseif(strcmpi(rescale_type{k}, 'output_rect'))
            ratefun = @(alpha, rateP, dat) alpha.*(exp(polyval(rateP, dat)) - exp(rateP(end))) + exp(rateP(end));
        elseif(strcmpi(rescale_type{k}, 'output'))
            ratefun = @(alpha, rateP, dat) alpha.*( exp(polyval(rateP, dat)) );
        end
        
        rates.fitfun{k} = ratefun;
        
        alpha_base = 1;
        nlogP = @(p, xt, xr) -sum(log(ratefun(alpha_base, p, xt) ) ) + sum(ratefun(alpha_base, p, xr)*deltaT);
        nlogP_resc = @(p, xt, xr) -sum(log(ratefun(p(end), p(1:end-1), xt) ) ) + sum(ratefun(p(end), p(1:end-1), xr)*deltaT);
        
        
        % fit High var. part first
        tic
        p0H = p0(tdatH, rdatH);
        nlogPH = @(p) nlogP(p, tdatH, rdatH);  % r(xH) = exp( a*x + b);
        [fP.fpH, ~] = fminunc(nlogPH, p0H, op); % p = [a b];
        toc
        
        %all - rescale Low Data
        tic
        nlogP_all = @(p) nlogP_resc(p, tdatL, rdatL) + nlogP(p(1:end-1), tdatH, rdatH);
        %         nlogP_all_fit = @(p) nlogP_all(p, tindsLfit, tindsHfit, rindsLfit, rindsHfit);
        
        % r(xH) = exp(a*xH + b); r(xL) = exp( a* (alphaL*xL) + b );
        [fP.rescL, LL.rescL(k)] = fminunc(nlogP_all, [fP.fpH, 1], op); % p = [a, b, alphaL];
        toc
        
        % fit nlogP to a quadratic function, and find the confidence interval for alphaL
        nLLfit = @(x, xdata) x(1)*(xdata-x(2)).^2 + x(3);
        alphas = linspace(1, 2.5, 30);
        ll = @(alpha) nlogP_all([fP.rescL(1:end-1), alpha]);
        for i=1:length(alphas), lla(i) = ll(alphas(i)); end
        aa = lsqcurvefit(nLLfit, [1, fP.rescL(1), -min(lla)], alphas, -lla);
        fP.s_alphaL = 1.645/sqrt(-2*aa(1));    % --- 90 percent confidence interval for alpha ( eqn. (2) in http://sites.stat.psu.edu/~abs12/stat504/Lecture/lec3_4up.pdf )
        
        return
        
        % fit low var. part first
        p0L = p0(tdatL, rdatL);
        nlogPL = @(p) nlogP(p, tdatL, rdatL);
        [fP.fpL, ~] = fminunc(nlogPL, p0L, op);
        
        % all - Rescale High Data
        
        % r(xH) = exp(a*(alphaH*xH) + b); r(xL) = exp( a*xL + b );
        nlogP_all = @(p, tlow, thigh, rlow, rhigh) nlogP_resc(p, tdatH(thigh), rdatH(rhigh)) + nlogP(p(1:end-1), tdatL(tlow), rdatL(rlow));
        nlogP_all_fit = @(p) nlogP_all(p, tindsLfit, tindsHfit, rindsLfit, rindsHfit);
        [fP.rescH, LL.rescH(k)] = fminunc(nlogP_all_fit, [fP.fpL, 1], op);

        
        %all no rescaling, sep. rates
        %     nlogP_noresc = @(p) nlogP(p(1:degree+1), tdata(tL), rdata(rL)) + nlogP(p(degree+2:end), tdata(tH), rdata(rH));
        %     [~, LL0] = fminunc(nlogP_noresc, [fpL fpH], op);
     
    end
   
end




