function [fitParams, logP, logP95] = fitRate(ratefun, x0, turndata, rundata, dt, nBins, histaxis)
%function [fitParams, logP, logP95] = fitRate(ratefun, x0, turndata, rundata, histaxis)
    

    
    nstim = size(turndata, 1);
    if(nstim>2)
        nstim = size(turndata, 2);
    end
    
    fitfun = @(x) -sum(log(max(ratefun(x, turndata), 1E-100))) + dt*sum(max(0,ratefun(x, rundata)));

    
    op = optimset('fmincon');
    op.Algorithm = 'interior-point';
    op.MaxFunEvals = 1e4;
    op.MaxIter = 1e4;
    %op.LargeScale = 'off';
    problem.solver = 'fmincon';
    problem.options = op;
    
    problem.x0 = x0;
    if(nstim==1) %1stim
        if(length(x0)==2)
            problem.lb = [0 -Inf];
            problem.ub = [10 Inf];
        else
            problem.lb = [-5 0 -Inf];
            problem.ub = [5 10 Inf];
        end
    elseif(nstim==2) %2stim
        if(length(x0)==2)
            problem.lb = [0 -Inf];
            problem.ub = [5 Inf];  
        elseif(length(x0)==3)
            problem.lb = [0 0 -Inf];
            problem.ub = [1 1 Inf];
        else
            problem.lb = [-1 -1 -1 0 0 -Inf];
            problem.ub = [1 1 1 3 3 Inf];
        end
    end

    problem.objective = fitfun;
    fitParams = fmincon(problem);
    logP = fitfun(fitParams);
    
%     ratefunOneParam = @(x, xdata) ratefun([x(1), fitParams(2)], xdata);
%     fitfunOneParam = @(x) -sum(log(max(ratefunOneParam(x, turndata), 1E-100))) + dt*sum(max(0,ratefunOneParam(x, rundata)));
%     
%     edge = 3.84;
%     xup = linspace(fitParams(1), fitParams(1)+1, 1e3);
%     xdown = linspace(fitParams(1)-3, fitParams(1), 1e3);
%     xx = [xdown xup];
%     
% 
%     for j=1:length(xx)
%         ratio(j) = -2*log(fitfunOneParam(xx(j)) / fitfunOneParam(fitParams(1)));
%     end
%     indup = find(ratio(end/2+1:end)<-edge, 1, 'first');
%     inddown = find(ratio(1:end/2)<-edge, 1, 'last');
%     edgeParams(1) = xup(indup);
%     edgeParams(2) = xdown(inddown);
%     logP95 = [fitfunOneParam(edgeParams(1)) fitfunOneParam(edgeParams(2))];
    
if(nstim==1)
    if (existsAndDefault('histaxis', []))
        ht = histcounts(turndata, binEdgesFromCenters(histaxis));
        hr = histcounts(rundata, binEdgesFromCenters(histaxis));
        r = ht./hr * 60/dt;
        eb = sqrt(ht)./hr * 60/dt;
        errorbar (histaxis, r, eb, 'ko');
        hold on;
        plot (histaxis, ratefun(fitParams, histaxis) *60);
%        plot(histaxis, 60*ratefunOneParam(edgeParams(1), histaxis), histaxis, 60*ratefunOneParam(edgeParams(2), histaxis));
        hold off
    end
elseif(nstim==2)
    if(size(turndata,2)==size(rundata,2))  
        adata = [turndata; rundata];
    else
        adata = [turndata, rundata];
    end
    
    if (existsAndDefault('histaxis', []))
        lx(:,1) = linspace(percentile(adata(:,1),0.01), percentile(adata(:,1),.99), nBins);
        lx(:,2) = linspace(percentile(adata(:,2),0.01), percentile(adata(:,2),.99), nBins);
        [lxx, lxy] = meshgrid(lx(:,1), lx(:,2));
        lxdata = [lxx(:) lxy(:)];
        h = makeIm(adata(:,1), adata(:,2), binEdgesFromCenters(lx(:,1)), binEdgesFromCenters(lx(:,2)));
        ht = makeIm(turndata(:,1), turndata(:,2), binEdgesFromCenters(lx(:,1)), binEdgesFromCenters(lx(:,2)));
        hist_all = reshape(h, size(lxx));
        hist_turn = reshape(ht, size(lxx));
        rate(:,1) = 60/dt * sum(hist_turn, 1)./sum(hist_all, 1);
        rate_eb(:,1) = 60/dt * sqrt(sum(hist_turn, 1))./sum(hist_all, 1);
        rate(:,2) = 60/dt * sum(hist_turn, 2)./sum(hist_all, 2);
        rate_eb(:,2) = 60/dt * sqrt(sum(hist_turn, 2))./sum(hist_all, 2);
        
        for k=1:2
            errorbar (lx(:,k), rate(:,k), rate_eb(:,k), 'marker', 'o', 'linestyle', '--');
            hold on;
            rateFit = reshape(ratefun(fitParams, lxdata), size(lxx));
            plot(lx(:,k), 60 * sum(rateFit.*hist_all, k)./sum(hist_all, k));
        end
    end
end
end