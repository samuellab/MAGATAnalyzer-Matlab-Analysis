function [ RF ] = plotScaledRates( struct, color, fit_test )


lograte = 1; %whether to plot rate functions on semi-log or not
deltaT=.05; % used to nornalize by animal-hours. this shouldn't be hard-coded, but it is for now

% Rescaled Rates

[adim,po] = blank8x10Figure();
pp = [adim.lx3, adim.h0 - adim.h, adim.w3, adim.h];
pos = {pp, pp, pp};
pos{2}(1) = adim.cx3;
pos{3}(1) = 2*pos{2}(1) - pos{1}(1);
pos{3}(2) = pos{1}(2);

pos2 = pos; % 2nd row of figure (3 columns each)
pos2{1}(1) = pos{1}(1);
pos2{1}(2) = adim.cx3 + 0.15;
pos2{2}(2) = pos2{1}(2);
pos2{3}(2) = pos2{1}(2);

pos3 = pos2;
pos3{1}(2) = 2*pos2{1}(2) - pos{1}(2);
pos3{2}(2) = pos3{1}(2);
pos3{3}(2) = pos3{1}(2);

pos4 = pos3;
pos4{1}(2) = 2*pos3{1}(2) - pos2{1}(2);
pos4{2}(2) = pos4{1}(2);
pos4{3}(2) = pos4{1}(2);

RedsHL = {1/255*[243 127 129], 1/255*[127 20 22]};
BluesHL = {1/255*[74 137 165], 1/255*[57 70 156]};

positions = {pos{1}, pos{2}, pos{3}, pos2{2}, pos3{2}};

if(strcmpi(color, 'red'))
    colors = RedsHL;
    col = [1 0 0];
    ymax = 12;
else
    colors = BluesHL;
    col = [0 0 1];
    ymax = 18;
end
barcolors = {[0 .5 0], [.2 .9 0], [0.5 1 0], col};

rates = struct.rates;

rPfit = [ mean(struct.adapt.quadR,2 ), mean(struct.adapt.linR,2 ), log(mean(struct.adapt.r0, 2)) ];
alphas = mean(struct.adapt.alpha, 2);

% plot fit-rate functions for each model

for k=1:length(struct.models)+1
    
    RF{k}.pos = positions{k}; %#ok<*AGROW>
    RF{k}.ax = axes('Position', positions{k});
    
    h(1) = errorbar (rates.lxL, 60*rates.rL, 60*rates.rL_eb, 'color', colors{1}, 'marker', 'o', 'linestyle', 'none');
    hold all;
    h(2) = errorbar (rates.lxH, 60*rates.rH, 60*rates.rH_eb, 'color', colors{2}, 'marker', 's','MarkerFaceColor', colors{2}, 'linestyle', 'none');
    
    
    % for k=1:length(struct.models)
    if(k>length(struct.models))
        ratefun = @(params,dat) exp(polyval(params, dat));
        rP_temp = reshape(cell2mat(struct.null.fP_noadapt), [length(struct.models) length(struct.null.nturns)]);
        rP = mean(rP_temp, 2);
        h(3) = plot (rates.lx_all, 60*ratefun(rP, rates.lx_all), 'color', colors{1}, 'linestyle', ':', 'LineWidth', 1.7);
        h(3) = plot (rates.lx_all, 60*ratefun(rP, rates.lx_all), 'color', colors{2}, 'linestyle', '-', 'LineWidth', 1);
    else
        ratefun = struct.ratefuns{k};
        alpha = alphas(k);
        rP = rPfit(k,:);
        
        h(3) = plot (rates.lxL, 60*ratefun(alpha, rP, rates.lxL), 'color', colors{1}, 'linestyle', ':', 'LineWidth', 1.7);
        h(4) = plot (rates.lxH, 60*ratefun(1, rP, rates.lxH), 'color', colors{2}, 'linestyle', '-', 'LineWidth', 1);
    end
    % end
    
    
    ax = gca;
    ax.Box = 'off';
    
    ax.YLabel.String = 'turn rate (min^{-1})';
    ax.YLabel.FontSize = 10;
    ax.XLabel.String = 'Filtered Stimulus Value';
    ax.XLabel.FontSize = 10;
    
    ax.XAxis.Color = col;
    
    RF{k}.handles = h;
    
    
    set([RF{k}.ax], 'YLim', [0.5 ymax]);
    if(lograte)
        set([RF{k}.ax], 'YScale', 'log');
        set([RF{k}.ax], 'YLim', [0.5 ymax]);
        set([RF{k}.ax], 'YTick', [0 1 2 5 10 15 20]);
        set([RF{k}.ax], 'YTickLabel', {'0', '1', '2', '5', '10', '15', '20'});
    end
    
end

k = k+1;


if(strcmpi(fit_test, 'fit'))
    
    RF{k}.pos = positions{k}; %#ok<*AGROW>
    RF{k}.ax = axes('Position', positions{k});

    LL_null = -struct.null.nLL_noadapt_fit;       %jwolk
%     LL_null2 = -struct.null.nLL_meanR;
    LLs_norm = -struct.adapt.nLL_fit - LL_null;
    
        deltaT=.05;
    norm = deltaT*struct.null.nruns/3600;
    LLs_norm = LLs_norm./norm;
    
    LLs = mean(LLs_norm,2);
    LLs_err = std(LLs_norm, 0, 2);%./sqrt(size(struct.adapt.nLL_fit, 2));
    
    for q=1:length(LLs)
        
        a = bar(q, LLs(q));
        a.BarWidth = .8;
        a.EdgeColor = col;
        a.FaceColor = col;
        a.FaceAlpha = 0.5;
        hold on;
        errorbar(q, LLs(q), LLs_err(q), 'color', col, 'LineWidth', 2, 'CapSize', 0);
    end
    ax = gca;
    RF{k}.ax = ax; 
    ax.XAxis.TickValues = [];
    
    
elseif(strcmpi(fit_test, 'test'))
    
    inds = struct.realinds_test;
    LL_null_fit = -struct.null.nLL_noadapt_fit(inds);
    LL_null_test = -struct.null.nLL_noadapt_test(inds);
    
    LL_fit = -struct.adapt.nLL_fit(:,inds) - LL_null_fit;
    LL_test = -struct.adapt.nLL_test(:,inds) - LL_null_test;
    
    norm = deltaT*struct.null.nruns(inds)/3600; % normalize by animal hours
    meanNorm = mean(norm);
    LL_test = LL_test./norm;
    LL_fit = LL_fit./norm;
    
    LL_fit = LL_fit.*meanNorm;
    LL_test = LL_test.*meanNorm;
    
    meanLL_fit = mean(LL_fit,2);
    errLL_fit= std(LL_fit, 0, 2)./sqrt(length(inds));
    meanLL_test = mean(LL_test,2);
    errLL_test = std(LL_test, 0, 2)./sqrt(length(inds));
    
    
    %log-likelihood for fit part, bar graph
%     RF{k}.pos = pos3{1}; %#ok<*AGROW>
%     RF{k}.ax = axes('Position', pos3{1});
%     for q=1:length(meanLL_fit)
%         
%         a = bar(q, meanLL_fit(q));
%         a.BarWidth = .8;
%         a.EdgeColor = col;
%         a.FaceColor = col;
%         a.FaceAlpha = 0.5;
%         hold on;
%         errorbar(q, meanLL_fit(q), errLL_fit(q), 'color', col, 'LineWidth', 2, 'CapSize', 0);
%     end
%     ax = gca;
%     RF{k}.ax = ax; 
%     ax.XAxis.TickValues = [];
    
    %log-likelihood for test part, bar graph
    RF{k}.pos = pos3{2}; %#ok<*AGROW>
    RF{k}.ax = axes('Position', pos3{2});
    for q=[1,2,3]
        
        a = bar(q, meanLL_test(q));
        a.BarWidth = .8;
        a.EdgeColor = col;
        a.FaceColor = col;
        a.FaceAlpha = 0.5;
        hold on;
        errorbar(q, meanLL_test(q), errLL_test(q), 'color', col, 'LineWidth', 2, 'CapSize', 0);
    end
    ax = gca;
    RF{k}.ax = ax;
    ax.XAxis.TickValues = [];
    ax.YAxis.Limits(1) = 0;
    
%     RF{k}.pos = pos4{1};
%     RF{k}.ax = axes('Position', pos4{1});
%     plot(LL_fit(3,:), LL_fit(1,:), 'k.', LL_fit(3,:), LL_fit(3,:), 'k--');
%     ax = gca;
%     axis equal;
    
    RF{k+1}.pos = pos4{2};
    RF{k+1}.ax = axes('Position', pos4{2});
    plot(LL_test(3,:), LL_test(1,:), 'k.', LL_test(3,:), LL_test(3,:), 'k--');
    ax = gca;
    
    RF{k+2}.pos = pos4{3};
    RF{k+2}.ax = axes('Position', pos4{3});
    histogram(LL_test(1,:)-LL_test(3,:), 15, 'FaceColor', 'k', 'FaceAlpha', 0.8);
    ax = gca;
%     axis equal;
%     for q=1:length(LLs)
%         
%         a = bar(q, LLs(q));
%         a.BarWidth = .8;
%         a.EdgeColor = col;
%         a.FaceColor = col;
%         a.FaceAlpha = 0.5;
%         hold on;
%         errorbar(q, LLs(q), LLs_err(q), 'color', col, 'LineWidth', 2, 'CapSize', 0);
%     end

%     a.Parent.YLim = [ round(min(LLs), -4), round(max(LLs), -4)];
    
end
    
    
    
    
    

