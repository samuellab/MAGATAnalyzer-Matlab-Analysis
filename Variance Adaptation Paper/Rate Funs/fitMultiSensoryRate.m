
function fitStruct = fitMultiSensoryRate( btdvar, nBins, pdegree, multiType, nExpTest)
% function fitStruct = fitMultiSensoryRate( btdvar, opstruct, multiType, nExpTest)

% very similar to fitRescaledRates.m, but this one is for 2-stimulus data:
% calculates fit parameters and log-likelihoods for 2-stim rate function


% nBins: number of bins the stimulus is divided into for rate-function plotting

% nExpTest is the number of experiments that are excluded from the fit

% multiType: {'ind', 'add', 'mult'} or just one of them, for r = r(xo)+r(xl), r=r(cos(th)*xo+sin(th)*xl) or r=r(xo)*r(xl) respectively

% pdegree: 1 or 2.
% note: if pdegree=1, add and mult are the same rate function

nexps = max(btdvar(1).fullensemble.expnum);
deltaT = median(diff(btdvar(1).fullensemble.eti));
theta0 = deg2rad(45);

op = optimoptions('fminunc');
op.Display = 'off';
op.Algorithm = 'quasi-newton';



for i=1:length(btdvar)
    varstruct = btdvar(i);
    td(:,i) = varstruct.turn.x_conv;
    rd(:,i) = varstruct.noturn.x_conv;
    teti = varstruct.turn.eti;
    reti = varstruct.noturn.eti;
end

fitStruct.Rates = turnRates2D(td, rd, deltaT, nBins);


if(pdegree==1)
    p0 = @(xt, xr) [mean(xt)/var(xr) log(length(xt)/(length(xr)*deltaT))];
else
    p0 = @(xt, xr) [-0.05 mean(xt)/var(xr) log(length(xt)/(length(xr)*deltaT))];
end




% fit all data to one rate function
% ratefun = @(rateP, dat) exp(polyval(rateP, dat));
% nlogP = @(p, xt, xr) -sum(log(ratefun(p, xt) ) ) + sum(ratefun(p, xr)*deltaT);
% fitfun = @(p) nlogP(p, tH_all, rH_all) + nlogP(p, tL_all, rL_all);
% fp0 = p0([tH_all;tL_all], [rH_all;rL_all]);
% [all.rP, all.nLL] = fminunc(fitfun, fp0, op);


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
    
    tnum = btdvar(1).turn.expnum;
    rnum = btdvar(1).noturn.expnum;
    
    ttest = tnum == expInds(qq,:);
    rtest = rnum == expInds(qq,:);
    fitExps = setdiff(1:nexps, expInds(qq,:));
    tfit = tnum == fitExps;
    rfit = rnum == fitExps;
    
    tdata1 = repmat(td(:,1), [1, size(tfit,2)]);
    tdata2 = repmat(td(:,2), [1, size(tfit,2)]);
    rdata1 = repmat(rd(:,1), [1, size(rfit,2)]);
    rdata2 = repmat(rd(:,2), [1, size(rfit,2)]);
    tdata_fit(:,1) = tdata1(tfit);
    tdata_fit(:,2) = tdata2(tfit);
    rdata_fit(:,1) = rdata1(rfit);
    rdata_fit(:,2) = rdata2(rfit);
    tdata1 = repmat(td(:,1), [1, size(ttest,2)]);
    tdata2 = repmat(td(:,2), [1, size(ttest,2)]);
    rdata1 = repmat(rd(:,1), [1, size(rtest,2)]);
    rdata2 = repmat(rd(:,2), [1, size(rtest,2)]);
    tdata_test(:,1) = tdata1(ttest);
    tdata_test(:,2) = tdata2(ttest);
    rdata_test(:,1) = rdata1(rtest);
    rdata_test(:,2) = rdata2(rtest);
    
    nturns_fit(i) = size(tdata_fit, 1);
    nruns_fit(i) = size(rdata_fit, 1);
    meanR_fit(i)= nturns_fit(i)/(nruns_fit(i)*deltaT);
    nturns_test(i) = size(tdata_test, 1);
    nruns_test(i) = size(rdata_test, 1);
    meanR_test(i)= nturns_test(i)/(nruns_test(i)*deltaT);
    
    % null model (mean rate)
    nLL_meanR_fit(i) = -nturns_fit(i)*log(meanR_fit(i)) + nruns_fit(i)*(meanR_fit(i)*deltaT);
    nLL_meanR_test(i) = -nturns_test(i)*log(meanR_test(i)) + nruns_test(i)*(meanR_test(i)*deltaT);
    
    if(pdegree==1)
        ratefun = @(rP, dat) exp( rP(1) * dat(:,1) + rP(2) * dat(:,2) + rP(3) );
        fitfun = ratefun;
        p01 = p0(tdata_fit(:,1), rdata_fit(:,1));
        p02 = p0(tdata_fit(:,2), rdata_fit(:,2));
        p0fit = [p01(1) p02];
        nlogP = @(p, xt, xr) -sum(log(ratefun(p, xt) ) ) + sum(ratefun(p, xr)*deltaT);
        
        nlogP_fit = @(p) nlogP(p, tdata_fit, rdata_fit);
        [fp, nLL_fit(i)] = fminunc(nlogP_fit, p0fit, op);
        nLL_test(i) = nlogP(fp, tdata_test, rdata_test) ;
        nLL_all(i) = nlogP(fp, [tdata_fit;tdata_test], [rdata_fit; rdata_test]);
        rateP{i} = fp;
        if(imag(nLL_test(i))==0)
            nLL_test_real = [nLL_test_real, nLL_test(:,i)];
        end
    elseif(pdegree==2)
        for k=1:length(multiType)
            if(strncmpi(multiType{k}, 'add', 3))
                ratefun = @(rP, dat) exp( rP(2) * (cos(rP(1))*dat(:,1) + sin(rP(1))*dat(:,2)).^2 + rP(3) * (cos(rP(1))*dat(:,1) + sin(rP(1))*dat(:,2)) + rP(4) );
                p01 = p0(tdata_fit(:,1), rdata_fit(:,1));
                p0fit = [theta0 0 0.5 p01(end)];
            elseif(strncmpi(multiType{k}, 'multiply', 3))
                ratefun = @(rP, dat) exp( rP(1)*dat(:,1).^2 + rP(2)*dat(:,1) + rP(3)*dat(:,2).^2 + rP(4)*dat(:,2) + rP(5) );
                p01 = p0(tdata_fit(:,1), rdata_fit(:,1));
                p02 = p0(tdata_fit(:,2), rdata_fit(:,2));
                p0fit = [p01(1:2), p02];
            elseif(strncmpi(multiType{k}, 'independent', 3))
                ratefun = @(rP, dat) exp( rP(1)*dat(:,1).^2 + rP(2)*dat(:,1) + rP(3) ) + exp( rP(4)*dat(:,2).^2 + rP(5)*dat(:,2) + rP(6) );
                p01 = p0(tdata_fit(:,1), rdata_fit(:,1));
                p02 = p0(tdata_fit(:,2), rdata_fit(:,2));
                p0fit = [p01, p02];
            end
            
            fitfun{k} = ratefun;
            
            nlogP = @(p, xt, xr) -sum(log(ratefun(p, xt) ) ) + sum(ratefun(p, xr)*deltaT);
            
            nlogP_fit = @(p) nlogP(p, tdata_fit, rdata_fit);
            [fp, nLL_fit(k,i)] = fminunc(nlogP_fit, p0fit, op);
            nLL_test(k, i) = nlogP(fp, tdata_test, rdata_test) ;
            nLL_all(k,i) = nlogP(fp, [tdata_fit;tdata_test], [rdata_fit; rdata_test]);
            rateP{k}{i} = fp;
            
            BIC_fit(k,i) = log(nturns_fit(i)+nruns_fit(i))*length(fp) + 2*nLL_fit(k,i);
            BIC_test(k,i) = log(nturns_test(i)+nruns_test(i))*length(fp) + 2*nLL_test(k,i);
            AIC(k,i) = 2*length(fp) + 2*nLL_fit(k,i);
            
        end
        if(any(imag(nLL_test(:,i))~=0)==0)
            nLL_test_real = [nLL_test_real, nLL_test(:,i)];
        end
    end
    

    
    if(size(nLL_test_real, 2)>=500)
        break
    end
    
    clear tdata_fit tdata_test rdata_fit rdata_test
    toc
end

for k=1:size(nLL_fit, 1)
    inds_test{k} = find(imag(nLL_test(k,:))==0);
    ninds(k) = length(inds_test{k});
end
if( any( ninds < size(expInds,1)) )
    m = find( ninds==min(ninds) );
    realinds_test = inds_test{m};
else
    realinds_test = inds_test{1};
end


fitStruct.models = multiType;
fitStruct.ratefuns = fitfun;
fitStruct.meanR_fit = meanR_fit;
fitStruct.meanR_test = meanR_test;
fitStruct.nturns_fit = nturns_fit;
fitStruct.nturns_test = nturns_test;
fitStruct.nruns_fit = nruns_fit;
fitStruct.nruns_test = nruns_test;
fitStruct.nLL_meanR_fit = nLL_meanR_fit;
fitStruct.nLL_meanR_test = nLL_meanR_test;
fitStruct.nLL_fit = nLL_fit;
fitStruct.nLL_test = nLL_test;
fitStruct.nLL_all = nLL_all;
fitStruct.BIC_fit = BIC_fit;
fitStruct.BIC_test = BIC_test;
fitStruct.AIC = AIC;
fitStruct.rP = rateP;
fitStruct.realinds_test = realinds_test;


end

function rf = turnRates2D(td, rd, deltaT, nBins)

lx(:,1) = linspace(-3*std(rd(:,1)), 3*std(rd(:,1)), nBins);
lx(:,2) = linspace(-3*std(rd(:,2)), 3*std(rd(:,2)), nBins);
[lxx, lxy] = meshgrid(lx(:,1), lx(:,2));
lxdata = [lxx(:) lxy(:)];

h = makeIm(rd(:,1), rd(:,2), binEdgesFromCenters(lx(:,1)), binEdgesFromCenters(lx(:,2)));
ht = makeIm(td(:,1), td(:,2), binEdgesFromCenters(lx(:,1)), binEdgesFromCenters(lx(:,2)));
rhist = reshape(h, size(lxx));
thist = reshape(ht, size(lxx));

rf.lx = lx;
rf.rhist = rhist;
rf.thist = thist;
rf.rate(:,1) = 1/deltaT * sum(thist, 1)./sum(rhist, 1);
rf.rate_eb(:,1) = 1/deltaT * sqrt(sum(thist, 1))./sum(rhist, 1);
rf.rate(:,2) = 1/deltaT * sum(thist, 2)./sum(rhist, 2);
rf.rate_eb(:,2) = 1/deltaT * sqrt(sum(thist, 2))./sum(rhist, 2);


end
