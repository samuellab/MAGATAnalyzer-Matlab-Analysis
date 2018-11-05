
function RescaleStruct = fitRescaledRates( btdvar, opstruct, degree, rescale_type, nExpTest )
% function RescaleStruct = fitRescaledRates( btdvar, opstruct, degree, rescale_type, nExpTest )

% jack-knife calculation of fit parameters and log-likelihood for fit to a rescaling model rate function

% degree = 1 or 2 for r=exp(ax+b) or r=exp(ax^2+bx+c) respectively
% rescale_type = {'input', 'output', 'output_rect'} or some subset of those

% nExpTest is the number of experiments from the total number that is excluded from the fits. 


op = optimoptions('fminunc');
op.Display = 'off';
op.Algorithm = 'quasi-newton';

nBins = opstruct.numLxBins;
deltaT = median(diff(btdvar.fullensemble.eti));

nexps = max(btdvar.fullensemble.expnum);

% find high and low rates from all the data
[tL_all, tH_all, rL_all, rH_all] = getLowHighData(btdvar, opstruct, 1:nexps);
rates.lxL = linspace(-3*std(rL_all), 3*std(rL_all), nBins);
rates.lxH = linspace(-3*std(rH_all), 3*std(rH_all), nBins);
rates.lx_all = linspace(-3*std([rL_all; rH_all]), 3*std([rL_all;rH_all]), nBins);
ntL = histcounts(tL_all, binEdgesFromCenters(rates.lxL));
nrL = histcounts(rL_all, binEdgesFromCenters(rates.lxL));
ntH = histcounts(tH_all, binEdgesFromCenters(rates.lxH));
nrH = histcounts(rH_all, binEdgesFromCenters(rates.lxH));
nt_all = histcounts([tL_all; tH_all], binEdgesFromCenters(rates.lx_all));
nr_all = histcounts([rL_all; rH_all], binEdgesFromCenters(rates.lx_all));
rates.rL = 1./deltaT *ntL./nrL;
rates.rL_eb = 1./deltaT * sqrt(ntL)./nrL;
rates.rH = 1./deltaT *ntH./nrH;
rates.rH_eb = 1./deltaT * sqrt(ntH)./nrH;
rates.rall = 1./deltaT*nt_all./nr_all;
rates.rall_eb = 1./deltaT * sqrt(nt_all)./nr_all;

RescaleStruct.rates = rates;

if(degree==1)
    p0 = @(xt, xr) [mean(xt)/var(xr) log(length(xt)/(length(xr)*deltaT))];
else
    p0 = @(xt, xr) [-0.05 mean(xt)/var(xr) log(length(xt)/(length(xr)*deltaT))];
end




% fit all data to one rate function
ratefun_all = @(rateP, dat) exp(polyval(rateP, dat));
nlogP = @(p, xt, xr) -sum(log(ratefun_all(p, xt) ) ) + sum(ratefun_all(p, xr)*deltaT);
fitfun_all = @(p) nlogP(p, tH_all, rH_all) + nlogP(p, tL_all, rL_all);
fp0 = p0([tH_all;tL_all], [rH_all;rL_all]);
[all.rP, all.nLL] = fminunc(fitfun_all, fp0, op);


expInds = nchoosek(1:nexps, nExpTest);

% if(size(expInds,1)>200)
%     n = 200;
% else
    n = size(expInds, 1);
% end

nLL_test_real = [];
array = 1:size(expInds,1);
for i=1:n
    tic
    qq = datasample(array, 1, 'Replace', false);
    array = setdiff(array, qq);
    [tL_test, tH_test, rL_test, rH_test] = getLowHighData(btdvar, opstruct, expInds(qq,:));
    [tL_fit, tH_fit, rL_fit, rH_fit] = getLowHighData(btdvar, opstruct, setdiff(1:nexps, expInds(qq,:)));
   
   nturns(i) = length([tL_fit;tH_fit]);
   nruns(i) = length([rL_fit;rH_fit]);
   meanR(i) = nturns(i)/(nruns(i)*deltaT);
   
   % null model (mean rate)
   nLL_meanR(i) = -nturns(i)*log(meanR(i)) + nruns(i)*(meanR(i)*deltaT);
   
   % no adaptation
   ratefun = @(rateP, dat) exp(polyval(rateP, dat));
   nlogP = @(p, xt, xr) -sum(log(ratefun(p, xt) ) ) + sum(ratefun(p, xr)*deltaT);
   nlogP_all = @(p, tlow, thigh, rlow, rhigh) nlogP(p, tlow, rlow) + nlogP(p, thigh, rhigh);
   nlogP_all_fit = @(p) nlogP_all(p, tL_fit, tH_fit, rL_fit, rH_fit);
   p0_noadapt = p0([tH_fit;tL_fit], [rH_fit;rL_fit]);
   
   [fp_noadapt{i}, nLL_noadapt_fit(i)] = fminunc(nlogP_all_fit, all.rP, op);
   nLL_noadapt_test(i) = nlogP_all(fp_noadapt{i}, tL_test, tH_test, rL_test, rH_test) ;
   
   for k=1:length(rescale_type)
        
        if(strcmpi(rescale_type{k}, 'input'))
            ratefun = @(alpha, rateP, dat) exp(polyval(rateP, alpha.*dat));
          elseif(strcmpi(rescale_type{k}, 'output_rect'))
            ratefun = @(alpha, rateP, dat) alpha.*(exp(polyval(rateP, dat)) - exp(rateP(end))) + exp(rateP(end));
        elseif(strcmpi(rescale_type{k}, 'output'))
            ratefun = @(alpha, rateP, dat) alpha.*( exp(polyval(rateP, dat)) );
        end
        
        fitfun{k} = ratefun;
        
        alpha_base = 1;
        nlogP = @(p, xt, xr) -sum(log(ratefun(alpha_base, p, xt) ) ) + sum(ratefun(alpha_base, p, xr)*deltaT);
        nlogP_resc = @(p, xt, xr) -sum(log(ratefun(p(end), p(1:end-1), xt) ) ) + sum(ratefun(p(end), p(1:end-1), xr)*deltaT);
        
%         % fit High var. part first
%         p0H = p0(tH_fit, rH_fit);
%         nlogPH = @(p) nlogP(p, tH_fit, rH_fit);
%         [fpH, ~] = fminunc(nlogPH, p0H, op);

        %all resc.
        nlogP_all = @(p, tlow, thigh, rlow, rhigh) nlogP_resc(p, tlow, rlow) + nlogP(p(1:end-1), thigh, rhigh);
        nlogP_all_fit = @(p) nlogP_all(p, tL_fit, tH_fit, rL_fit, rH_fit);
        
        [fp, nLL_fit(k,i)] = fminunc(nlogP_all_fit, [all.rP, 1], op);
        nLL_test(k, i) = nlogP_all(fp, tL_test, tH_test, rL_test, rH_test) ;
        nLL_all(k,i) = nlogP_all(fp, [tL_fit;tL_test], [tH_fit;tH_test], [rL_fit;rL_test], [rH_fit;rH_test]);
        if(degree==2)
            rP_quad(k,i) = fp(1);
        end
        rP_lin(k,i) = fp(end-2);
        rP_r0(k,i) = fp(end-1);
        alpha_resc(k,i) = fp(end);
    
    
   end
   if(any(imag(nLL_test(:,i))~=0)==0)
       nLL_test_real = [nLL_test_real, nLL_test(:,i)];
   end
   
%    if(size(nLL_test_real, 2)>=200)
%        break
%    end
   
   toc
end





for k=1:length(rescale_type)
    inds_test{k} = find(imag(nLL_test(k,:))==0);
    ninds(k) = length(inds_test{k});
end
if( any( ninds < size(expInds,1)) )
    m = find( ninds==min(ninds) );
    realinds_test = inds_test{m};
else
    realinds_test = inds_test{1};
end

RescaleStruct.models = rescale_type;
RescaleStruct.ratefuns = fitfun;
RescaleStruct.imaginds_test = inds_test;
RescaleStruct.realinds_test = realinds_test;
RescaleStruct.null.nturns = nturns;
RescaleStruct.null.nruns = nruns;
RescaleStruct.null.meanR = meanR;
RescaleStruct.null.nLL_meanR = nLL_meanR;
RescaleStruct.null.fP_noadapt = fp_noadapt;
RescaleStruct.null.nLL_noadapt_fit = nLL_noadapt_fit;
RescaleStruct.null.nLL_noadapt_test = nLL_noadapt_test;
RescaleStruct.adapt.nLL_fit = nLL_fit;
RescaleStruct.adapt.nLL_test = nLL_test;
RescaleStruct.adapt.nLL_all = nLL_all;
RescaleStruct.adapt.alpha = alpha_resc;
RescaleStruct.adapt.linR = rP_lin;
if(degree==2)
   RescaleStruct.adapt.quadR = rP_quad;
end
RescaleStruct.adapt.r0 = exp(rP_r0);


end

% 
% for i=1:n
%     tic
%         
%     [rHfit, rindsHfit] = datasample(rdatH, round(jack_knife*length(rdatH)), 'Replace', false);
%     [tHfit, tindsHfit] = datasample(tdatH, round(jack_knife*length(tdatH)), 'Replace', false);
%     [rLfit, rindsLfit] = datasample(rdatL, round(jack_knife*length(rdatL)), 'Replace', false);
%     [tLfit, tindsLfit] = datasample(tdatL, round(jack_knife*length(tdatL)), 'Replace', false);
%     if(jack_knife==1)
%         rindsHtest = rindsHfit;
%         rindsLtest = rindsLfit;
%         tindsHtest = tindsHfit;
%         tindsLtest = tindsLfit;
%     else
%         rindsHtest = setdiff(1:length(rdatH),rindsHfit);
%         tindsHtest = setdiff(1:length(tdatH),tindsHfit);
%         rindsLtest = setdiff(1:length(rdatL),rindsLfit);
%         tindsLtest = setdiff(1:length(tdatL),tindsLfit);
%     end
%     
%     for k=1:length(rescale_type)
%         
%         if(strcmpi(rescale_type{k}, 'input'))
%             ratefun = @(alpha, rateP, dat) exp(polyval(rateP, alpha.*dat));
%         elseif(strcmpi(rescale_type{k}, 'output_rect'))
%             ratefun = @(alpha, rateP, dat) alpha.*(exp(polyval(rateP, dat)) - exp(rateP(end))) + exp(rateP(end));
%         elseif(strcmpi(rescale_type{k}, 'output'))
%             ratefun = @(alpha, rateP, dat) alpha.*( exp(polyval(rateP, dat)) );
%         end
%         
%         rates.fitfun{k} = ratefun;
%         
%         alpha_base = 1;
%         nlogP = @(p, xt, xr) -sum(log(ratefun(alpha_base, p, xt) ) ) + sum(ratefun(alpha_base, p, xr)*deltaT);
%         nlogP_resc = @(p, xt, xr) -sum(log(ratefun(p(end), p(1:end-1), xt) ) ) + sum(ratefun(p(end), p(1:end-1), xr)*deltaT);
%         
%         
%         % fit High var. part first
%             
%         p0H = p0(tHfit, rHfit);
%         nlogPH = @(p) nlogP(p, tHfit, rHfit);
%         [fpH, ~] = fminunc(nlogPH, p0H, op);
%         
%         
%         %all resc.
%         
%         nlogP_all = @(p, tlow, thigh, rlow, rhigh) nlogP_resc(p, tdatL(tlow), rdatL(rlow)) + nlogP(p(1:end-1), tdatH(thigh), rdatH(rhigh));
%         
%         nlogP_all_fit = @(p) nlogP_all(p, tindsLfit, tindsHfit, rindsLfit, rindsHfit);
%         [fp, LL_fit(k,i)] = fminunc(nlogP_all_fit, [fpH, 1], op);
%         LL_test(k,i) = nlogP_all(fp, tindsLtest, tindsHtest, rindsLtest, rindsHtest) ;
%         LL_all(k,i) = nlogP_all(fp, [tindsLtest tindsLfit], [tindsHtest tindsHfit], [rindsLtest rindsLfit], [rindsHtest rindsHfit]);
%         
%         if(degree==2)
%             rP_quad(k,i) = fp(1);
%         end
%         rP_lin(k,i) = fp(end-2);
%         rP_r0(k,i) = fp(end-1);
%         alpha_resc(k,i) = fp(end);
%         
%         
%         %         nlogP_all_test = @(p) nlogP_all(p, tindsLtest, tindsHtest, rindsLtest, rindsHtest);
%         %         [fp_test{k}{i}, LL_test.fit(k,i)] = fminunc(nlogP_all_test, [fpH, 1], op);
% %         LL_test.test(k,i) = nlogP_all(fp_test{k}{i}, tindsLfit, tindsHfit, rindsLfit, rindsHfit) ;
%         
%         % fit low var. part first
%         % p0L = p0(tdatL(tindsL), rdatL(rindsL));
%         % nlogPL = @(p) nlogP(p, tdatL(tindsL), rdatL(rindsL));
%         % [fpL, llL] = fminunc(nlogPL, p0L, op);
%         
%         % nlogP_all = @(p) nlogP_resc(p, tdata(tH), rdata(rH)) + nlogP(p(1:end-1), tdata(tL), rdata(rL));
%         % [fp, LL.LtoH] = fminunc(nlogP_all, [fpL, 1], op);
%         %
%         % rP.LtoH = fp(1:end-1);
%         % alpha.LtoH = fp(end);
%         
%         %all no rescaling, sep. rates
%         %     nlogP_noresc = @(p) nlogP(p(1:degree+1), tdata(tL), rdata(rL)) + nlogP(p(degree+2:end), tdata(tH), rdata(rH));
%         %     [~, LL0] = fminunc(nlogP_noresc, [fpL fpH], op);
%      
%     end
%    
%     
%     toc
% end
% 
% for k=1:length(rescale_type)
%     
%     inds{k} = find(imag(LL_test(k,:))==0);
%     ninds(k) = length(inds{k});
% end
% if( any( ninds < n*length(rescale_type)) )
%     m = find( ninds==min(ninds) );
%     realinds = inds{m};
% end
% 
% JK.nLL_fit = LL_fit(:,realinds);
% JK.nLL_test = LL_test(:,realinds);
% JK.nLL_all = LL_all(:,realinds);
% JK.alpha = alpha_resc(:,realinds);
% JK.rP.lin = rP_lin(:,realinds);
% if(degree==2)
%     JK.rP.quad = rP_quad(:,realinds);
% end
% JK.rP.r0 = rP_r0(:,realinds);
% end


function [tL, tH, rL, rH] = getLowHighData(btdvar, opstruct, expInds)


tnum = btdvar.turn.expnum;
rnum = btdvar.noturn.expnum;

tvalid = tnum == expInds;
rvalid = rnum == expInds;

tdata = repmat(btdvar.turn.x_conv, [1, size(tvalid,2)]);
rdata = repmat(btdvar.noturn.x_conv, [1, size(rvalid,2)]);
tdata = tdata(tvalid);
rdata = rdata(rvalid);
    
tton = repmat(btdvar.turn.ton, [1, size(tvalid,2)]); tton = tton(tvalid);
ttoff = repmat(btdvar.turn.toff, [1, size(tvalid,2)]); ttoff = ttoff(tvalid);
teti = repmat(btdvar.turn.eti, [1, size(tvalid,2)]); teti = teti(tvalid);
rton = repmat(btdvar.noturn.ton, [1, size(rvalid,2)]); rton = rton(rvalid);
rtoff = repmat(btdvar.noturn.toff, [1, size(rvalid,2)]); rtoff = rtoff(rvalid);
reti = repmat(btdvar.noturn.eti, [1, size(rvalid,2)]); reti = reti(rvalid);


tH = tton > opstruct.adaptationTime & tton < ttoff & teti >= min(opstruct.timeRange) & teti < max(opstruct.timeRange);
tL = ttoff > opstruct.adaptationTime & ttoff < tton & teti >= min(opstruct.timeRange) & teti < max(opstruct.timeRange);
tU = tton > 0 & tton < opstruct.adaptationTime & teti >= min(opstruct.timeRange) & teti < max(opstruct.timeRange);
tD = ttoff > 0  & ttoff < opstruct.adaptationTime & teti >= min(opstruct.timeRange) & teti < max(opstruct.timeRange);
rH = rton > opstruct.adaptationTime & rton < rtoff & reti >= min(opstruct.timeRange) & reti < max(opstruct.timeRange);
rL = rtoff > opstruct.adaptationTime & rtoff < rton & reti >= min(opstruct.timeRange) & reti < max(opstruct.timeRange);
rU = rton > 0 & rton < opstruct.adaptationTime & reti >= min(opstruct.timeRange) & reti < max(opstruct.timeRange);
rD = rtoff > 0  & rtoff < opstruct.adaptationTime & reti >= min(opstruct.timeRange) & reti < max(opstruct.timeRange);


tH = tdata(tH);
tL = tdata(tL);
rH = rdata(rH);
rL = rdata(rL);
end



