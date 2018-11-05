function ah = VarAdaptPlots (btdstruct, names, colors, pos, exptype, deltaT, estimator, simtype, dt_opt, tau_opt, params_opt, StimNum, pdegree, multiType)
% function ah = VarAdaptPlots (btdstruct, names, colors, pos, exptype, deltaT, estimator, simtype, dt_opt, tau_opt, params_opt, StimNum, pdegree, multiType)

% see VarAdaptFigures for examples on how to use this script to plot stuff


switch exptype
    
    
    case '2stimRF'
        
        for j=1:length(btdstruct)
            
            if(isfield(btdstruct(j), 'RatesU'))
                Rates = btdstruct(j).RatesU;
            elseif(isfield(btdstruct(j), 'Rates'))
                Rates = btdstruct(j).Rates;
            else
                Rates = btdstruct(j);
            end
            
            
            
            lx_low = Rates.rfLowVar(1).lx;
            lx_high = Rates.rfHighVar(1).lx;
            rateX_low = 60*Rates.rfLowVar(1).rate;
            rateX_low_eb = 60*Rates.rfLowVar(1).rate_eb;
            rateX_high = 60*Rates.rfHighVar(1).rate;
            rateX_high_eb = 60*Rates.rfHighVar(1).rate_eb;
            
            ly_low = Rates.rfLowVar(2).lx;
            ly_high = Rates.rfHighVar(2).lx;
            rateY_low = 60*Rates.rfLowVar(2).rate;
            rateY_low_eb = 60*Rates.rfLowVar(2).rate_eb;
            rateY_high = 60*Rates.rfHighVar(2).rate;
            rateY_high_eb = 60*Rates.rfHighVar(2).rate_eb;
            
            if(pdegree==2)
                rateX_fit_low = 60*exp(polyval(Rates.rfLowVar(1).rateFitQuad, lx_low));
                rateX_fit_high = 60*exp(polyval(Rates.rfHighVar(1).rateFitQuad, lx_high));
                rateY_fit_low = 60*exp(polyval(Rates.rfLowVar(2).rateFitQuad, ly_low));
                rateY_fit_high = 60*exp(polyval(Rates.rfHighVar(2).rateFitQuad, ly_high));
            else
                rateX_fit_low = 60*exp(polyval(Rates.rfLowVar(1).rateFitLin, lx_low));
                rateX_fit_high = 60*exp(polyval(Rates.rfHighVar(1).rateFitLin, lx_high));
                rateY_fit_low = 60*exp(polyval(Rates.rfLowVar(2).rateFitLin, ly_low));
                rateY_fit_high = 60*exp(polyval(Rates.rfHighVar(2).rateFitLin, ly_high));
            end
            
            
            ah{1}(j).pos = pos{j}{1}; %#ok<*AGROW>
            ah{1}(j).ax = axes('Position', pos{j}{1});
            ah{1}(j).name = names{j};
            
            
            h(1) = errorbar (lx_low, rateX_low, rateX_low_eb, 'color', colors{1}{1}, 'marker', 'o', 'linestyle', 'none');
            hold all;
            h(2) = plot (lx_low, rateX_fit_low, 'color', colors{1}{1}, 'linestyle', '--', 'LineWidth', 2);
            
            h(3) = errorbar (lx_high, rateX_high, rateX_high_eb, 'color', colors{1}{2}, 'marker', 's', 'MarkerFaceColor', colors{1}{2}, 'linestyle', 'none');
            h(4) = plot (lx_high, rateX_fit_high, 'color', colors{1}{2}, 'linestyle', '-', 'LineWidth', 2);
            hold off;
            
            ah{1}(j).handles = h;
            title(names{j});
            ax = gca;
            ax.Box = 'off';
            %             set(ax, 'YScale', 'log');
            ax.YLim = [0 15];
            ax.YTick = [0 5 10 15 20];
            ax.YLabel.String = 'turn rate (min^{-1})';
            ax.YLabel.FontSize = 10;
            ax.XLabel.String = 'Filtered Stimulus Value';
            ax.XLabel.FontSize = 10;
            ax.XAxis.Color = 'r';
            
            ah{2}(j).pos = pos{j}{2}; %#ok<*AGROW>
            ah{2}(j).ax = axes('Position', pos{j}{2});
            ah{2}(j).name = names{j};
            
            hold on
            
            h2(1) = errorbar (ly_low, rateY_low, rateY_low_eb, 'color', colors{2}{1}, 'marker', 'o', 'linestyle', 'none');
            hold all;
            h2(2) = plot (ly_low, rateY_fit_low, 'color', colors{2}{1}, 'linestyle', '--', 'LineWidth', 2);
            
            h2(3) = errorbar (ly_high, rateY_high, rateY_high_eb, 'color', colors{2}{2}, 'marker', 's', 'MarkerFaceColor', colors{2}{2}, 'linestyle', 'none');
            h2(4) = plot (ly_high, rateY_fit_high, 'color', colors{2}{2}, 'linestyle', '-', 'LineWidth', 2);
            
            ah{2}(j).handles = h2;
            title(names{j});
            ax = gca;
            ax.Box = 'off';
            %             set(ax, 'YScale', 'log');
            ax.YLim = [0 15];
            ax.YTick = [0 5 10 15 20];
            ax.YLabel.String = 'turn rate (min^{-1})';
            ax.YLabel.FontSize = 10;
            ax.XLabel.String = 'Filtered Stimulus Value';
            ax.XLabel.FontSize = 10;
            ax.XAxis.Color = 'b';
            
            %         set(ax, 'YScale', 'log');
            
            
            
            
        end
        
    case 'ScaleFactorVTime'
        
        op = optimoptions('fminunc');
        op.Algorithm = 'quasi-newton';
        problem.options = op;
        problem.solver = 'fminunc';
        
        for j=1:length(btdstruct)
            
            ah(j).pos = pos{j}; %#ok<*AGROW>
            ah(j).ax = axes('Position', pos{j});
            axis square
            ah(j).name = names{j};
            
            clear tSim aSim sSim vSim tEst aEst vEst
            
            deltaT = median(diff(btdstruct(j).var(1).fullensemble.eti));

            if(strcmpi(estimator, 'PPF'))
                if(isfield(btdstruct, 'var_uv'))
                    if(StimNum == 1)
                        alphastruct = btdstruct(j).alphavseti_U;
                    else
                        alphastruct = btdstruct(j).alphavseti_UV;
                    end
                else
                    alphastruct = btdstruct(j).alphavseti;
                end
            elseif(strcmpi(estimator, 'Bayes'))
                if(isfield(btdstruct, 'var_uv'))
                    if(StimNum == 1)
                        alphastruct = btdstruct(j).BayesAlphavsEti_U;
                    else
                        alphastruct = btdstruct(j).BayesAlphavsEti_UV;
                    end
                else
                    alphastruct = btdstruct(j).BayesAlphavsEti_OL;
                end
                
            end
            
            tx_ton = alphastruct.tx_ton;
            alpha_ton = alphastruct.alpha_ton;
            valpha_ton = sqrt(alphastruct.valpha_ton);
            if(size(alpha_ton,1)>size(alpha_ton,2))
                alpha_ton = alpha_ton';
            end
            nstim = size(alpha_ton, 1);
            
            tshift = 3;
            
            if(nstim==2)
                period = btdstruct(j).var_uv(1).period;
            else
                period = btdstruct(j).var.period;
            end
            
            dt = median(diff(tx_ton));
            
            for k=1:size(alpha_ton, 1)
                alpha_shift(k,:) = circshift(alpha_ton(k,:), round(length(tx_ton)/2));
                if(nstim==2)
                    valpha_shift(k,:) = circshift(valpha_ton(k, k,:), round(length(tx_ton)/2));
                else
                    valpha_shift(k,:) = circshift(valpha_ton(k,:), round(length(tx_ton)/2));
                end
                aEst(k,:) = [alpha_shift(k, round(end-tshift/dt):end), alpha_shift(k,:)];
                vEst(k,:) = [valpha_shift(k, round(end-tshift/dt):end), valpha_shift(k,:)];
            end
            tx_shift = -fliplr(tx_ton);
            %             tEst = [-tshift:dt:tx_ton(1), tx_ton];
            tEst = [tx_shift(round(end-tshift/dt):end), tx_ton];
            if(max(tEst)>50)
                ind = find(tEst>50, 1);
                tEst = tEst(1:ind);
                aEst = aEst(k, 1:ind);
                vEst = vEst(k, 1:ind);
            end
            
            if(strcmpi(simtype, 'LowHigh'))
                                
                sim = btdstruct(j).LowHighSim;
                if(strcmpi(estimator, 'PPF'))
                    txSim = sim.alphavseti{end}.tx_ton;
                    alphaSim =  sim.alphavseti{end}.alpha_ton;
                    valphaSim = sqrt(sim.alphavseti{end}.valpha_ton);
                elseif(strcmpi(estimator, 'Bayes'))
                    txSim = sim.BayesAlphavsEti{end}.tx_ton;
                    alphaSim =  sim.BayesAlphavsEti{end}.alpha_ton;
                    valphaSim = sqrt(sim.BayesAlphavsEti{end}.valpha_ton);
                end
                dt = median(diff(txSim));
                tx_shift = circshift(txSim, round(length(txSim/2)));
                alphaSim_shift = circshift(alphaSim, round(length(txSim)/2));
                valphaSim_shift = circshift(valphaSim, round(length(txSim)/2));
                
                tSim = [-(txSim(1:round(1+tshift/dt))), txSim];
                aSim = [alphaSim_shift(round(end-tshift/dt):end), alphaSim_shift];
                aSim = aSim./mean(aSim);
                vSim = [valphaSim_shift(round(end-tshift/dt):end), valphaSim_shift];
                
            elseif(strcmpi(simtype, 'Optimal'))
                
                isconvolved = 0;
                
                tton = btdstruct(j).var.turn.ton;
                rton = btdstruct(j).var.noturn.ton;
                teti = btdstruct(j).var.turn.eti;
                reti = btdstruct(j).var.noturn.eti;
                xt = btdstruct(j).var.turn.x_conv;
                xr = btdstruct(j).var.noturn.x_conv;
                
%                 gainFun = @(x, xdata) 1./(sqrt(xdata.^2 + x(1)))./mean(1./(sqrt(xdata.^2 + x(1))));
                gainFun = @(x, xdata) 1./(sqrt(xdata.^2 + x^2));
                VarEst = @(dt, tau) AlphaEstimate(btdstruct(j), btdstruct(j).varops, 'eti', isconvolved, dt/tau, dt, 60:dt:1200, btdstruct(j).var.period, btdstruct(j).var.tshift);
                
                
                for i=1:length(dt_opt)
                    tic
                    sim{i} = VarEst(dt_opt(i), tau_opt(i));
                    toc
                end

                for k=1:length(sim)
                    
                    tx_ton = sim{k}.tx_ton;
                    txSim{k} = sim{k}.tx_ton;
                    tx = sim{k}.tx;
                    dt = median(diff(txSim{k}));
                    sigma = sim{k}.sigma;
                    sigma_ton = sim{k}.sigma_ton;
                    ssigma_ton = sqrt(sim{k}.vsigma_ton);
                    
                    fpar = params_opt{j};
                    
                    data = alphastruct.scaledValues;
                    
                    alphaup.alpha_ton = gainFun(fpar, sigma_ton+ssigma_ton);
                    alphaup.tx_ton = tx_ton;
                    alphaup = normAlpha(alphaup, data, 'ton', []);
                    
                    alpha.tx_ton = tx_ton;
                    alpha.alpha_ton = gainFun(fpar, sigma_ton);
                    alpha = normAlpha(alpha, data, 'ton', []);
                    
                    alphaSim{k} = alpha.alpha_ton;
                    salphaSim{k} = alphaup.alpha_ton - alpha.alpha_ton;
                    
                    alphaSim_shift{k} = circshift(alphaSim{k}, round(length(txSim{k})/2));
                    salphaSim_shift{k} = circshift(salphaSim{k}, round(length(txSim{k})/2));
                    aSim{k} = [alphaSim_shift{k}(round(end-tshift/dt):end), alphaSim_shift{k}];
                    npoints = length(aSim{k}) - length(alphaSim{k});
                    sSim{k} = [salphaSim_shift{k}(round(end-tshift/dt):end), salphaSim_shift{k}];
                    tSim{k} = [txSim{k}(end-npoints+1:end)-period, txSim{k}];
                end
                
            end
            
            
            
            if(isempty(estimator))
                
                temp = shadedErrorPlot({tEst tEst}, mat2cell(aEst, [1 1], length(aEst)), mat2cell(vEst, [1 1], length(vEst)));
                h(1) = temp(1);
                h(2) = temp(2);
                h(3) = temp(3);
                h(4) = temp(4);
                pp = h(1).Parent;
                pp.Children(1).Color = colors{j}{1}{1};
                pp.Children(2).Color = colors{j}{1}{2};
                pp.Children(1).LineWidth = 1.5;
                pp.Children(2).LineWidth = 1.5;
                pp.Children(3).FaceColor = colors{j}{1}{2};
                pp.Children(4).FaceColor = colors{j}{2}{2};
                pp.Children(3).FaceAlpha = .5;
                pp.Children(4).FaceAlpha = 1;
                
                %                 uistack(pp.Children(2),'top');
                
            elseif(strcmpi(simtype, 'LowHigh'))
                
                temp = shadedErrorPlot({tEst tSim}, {aEst aSim}, {vEst vSim});
                h(1) = temp(1);
                h(2) = temp(2);
                h(3) = temp(3);
                h(4) = temp(4);
                pp = h(1).Parent;
                pp.Children(1).Color = colors{j}{1}{1};
                pp.Children(2).Color = [0.1 0.1 0.1 .3];
                pp.Children(1).LineWidth = 1.5;
                pp.Children(2).LineWidth = 1.2;
                pp.Children(3).FaceColor = colors{j}{1}{2};
                pp.Children(4).FaceColor = [0.6 0.6 0.6];
                pp.Children(3).FaceAlpha = 1;
                pp.Children(4).FaceAlpha = .3;
                
            elseif(strcmpi(simtype, 'Optimal'))
                
                
                if(strcmpi(btdstruct(j).var.gqname, 'led1ValDiff'))
                    calpha = autumn;
                else
                    calpha = winter;
                end
                
                %                 calpha = winter;
                
                N = length(tSim);
                tSim{end+1} = tEst;
                aSim{end+1} = aEst;
                sSim{end+1} = vEst;
                temp = shadedErrorPlot(tSim, aSim, sSim);
                
                h(1) = temp(1);
                h(2) = temp(2);
                h(3) = temp(3);
                h(4) = temp(4);
                pp = h(1).Parent;
                
                nn = floor(length(calpha)/N);
                for k=1:N
                    pp.Children(k).Color = calpha(nn*k,:);
                    cols(k,:) = calpha(nn*k,:);
                    pp.Children(k).LineWidth = 1.5;
                    pp.Children(k+1+N).FaceColor = calpha(nn*k,:);
                    pp.Children(k+1+N).FaceAlpha = .2;
                end
                %pp.Children(N+1).Color = colors{j}{1}{1};
                pp.Children(N+1).Color = 'k';
                pp.Children(N+1).LineWidth = 1.5;
                %pp.Children(end).FaceColor = colors{j}{1}{2};
                pp.Children(end).FaceColor = [.6 .6 .6];
                pp.Children(end).FaceAlpha = .6;
                
                ax = gca;
                if(std(tau_opt)>std(dt_opt))
                    paraxis = tau_opt;
                else
                    paraxis = dt_opt;
                end
                ax.CLim = [min(paraxis) max(paraxis)];
                a{j} = colorbar;
                a{j}.Parent.Colormap = cols;
                lims = a{j}.Limits;
                int = diff(lims)/N;
                for q=1:N, xx(q) = lims(1) + q*int - int/2; end
                a{j}.Ticks = xx;
                a{j}.TickLabels = strsplit(num2str(paraxis));
                a{j}.TickLength = 0;
                
            else
                
                if(nstim==2)
                    h = shadedErrorPlot({tEst tEst}, {aEst(1,:) aEst(2,:)},{vEst(1,:) vEst(2,:)}, {vEst(1,:) vEst(2,:)}, {colors{j}{1}{1}, colors{j}{1}{2}});
                    pp = h.Parent;
                    pp.Children(1).LineWidth = 1.5;
                    pp.Children(2).LineWidth = 1.5;
                else
                    temp = shadedErrorPlot(tEst, aEst, vEst);
                    h(1) = temp(1);
                    h(2) = temp(2);
                    pp = h(1).Parent;
                    pp.Children(1).Color = colors{j}{1}{1};
                    pp.Children(1).LineWidth = 1.5;
                    pp.Children(2).FaceColor = colors{j}{1}{2};
                end
                
                
            end
            %
            
            ax = gca;
            ax.YLim = [0.6 1.6];
            ax.YAxis(1).Limits = ax.YLim;
            hold on
            h(end+1) = plot([0 0], [ax.YLim(1) ax.YLim(2)], 'k--', 'linewidth', 1.5);
            h(end+1) = plot([period/2 period/2], [ax.YLim(1) ax.YLim(2)], 'k--', 'linewidth', 1.5);
            hold off
            
            ax.YAxis(1).Color = colors{j}{1}{1};
            ax.YLabel.String = '\alpha';
            ax.YLabel.FontSize = 10;
            ax.XLabel.String = 'time in cycle (s)';
            ax.XLabel.FontSize = 10;
            
            
            ARatio = .66;
            
            ax = h(1).Parent;
            ax.DataAspectRatio(1) = ARatio*ax.DataAspectRatio(1);
            ax.Box = 'off';
            
            
            
            
        end
        
    case '1stimRF'
        
        SF = 60;
        
        
        for j = 1:length(btdstruct)
            
            
            if( strncmpi(btdstruct(j).var(StimNum).gqname, 'led2', 4))
                stim = 'blue';
            elseif( strncmpi(btdstruct(j).var(StimNum).gqname, 'led1', 4))
                stim = 'red';
            end
            
            ah(j).pos = pos{j}; %#ok<*AGROW>
            ah(j).ax = axes('Position', pos{j});
            ah(j).name = names{j};
            
            
            if(isfield(btdstruct(j), 'RatesU'))
                Rates = btdstruct(j).RatesU;
                stim = 'purple';
            elseif(isfield(btdstruct(j), 'Rates'))
                Rates = btdstruct(j).Rates;
            else
                Rates = btdstruct(j);
            end
            
            if(nargin>5)
                k = StimNum;
            else
                k=1;
            end
            
            lx_low = Rates.rfLowVar(k).lx;
            lx_high = Rates.rfHighVar(k).lx;
            rate_low = SF*Rates.rfLowVar(k).rate;
            rate_eb_low = SF*Rates.rfLowVar(k).rate_eb;
            rate_high = SF*Rates.rfHighVar(k).rate;
            rate_eb_high = SF*Rates.rfHighVar(k).rate_eb;
            
            if(pdegree==1)
                
                fitRate_low = SF*exp(polyval(Rates.rfLowVar(k).rateFitLin, lx_low));
                fitRate_high = SF*exp(polyval(Rates.rfHighVar(k).rateFitLin, lx_high));
            elseif(pdegree==2)
                fitRate_low = SF*exp(polyval(Rates.rfLowVar(k).rateFitQuad, lx_low));
                fitRate_high = SF*exp(polyval(Rates.rfHighVar(k).rateFitQuad, lx_high));
            end
            
            h(1) = errorbar (lx_low, rate_low, rate_eb_low, 'color', colors{j}{1}, 'marker', 'o', 'linestyle', 'none'); %'MarkerSize', 5, 'CapSize', 5
            hold all;
            h(2) = errorbar (lx_high, rate_high, rate_eb_high, 'color', colors{j}{2}, 'marker', 's', 'MarkerFaceColor', colors{j}{2}, 'linestyle', 'none');
            
            if(pdegree>0)
                h(3) = plot (lx_low, fitRate_low, 'color', colors{j}{1}, 'linestyle', ':', 'LineWidth', 2);
                h(4) = plot (lx_high, fitRate_high, 'color', colors{j}{2}, 'linestyle', '-', 'LineWidth', 2);
            end
            
            title(names{j});
            
            ax = gca;
            ax.Box = 'off';
            %             axis square
            
            ax.YLim = [0 15];
            % ax.YTick = [0 5 10];
            ax.YLabel.String = 'turn rate (min^{-1})';
            ax.YLabel.FontSize = 10;
            ax.XTick = [-10 -5 0 5 10];
            ax.XLabel.String = 'Filtered Stimulus Value';
            ax.XLabel.FontSize = 10;
            
            if(strcmpi(stim, 'red'))
                ax.XAxis.Color = 'r';
            elseif(strcmpi(stim, 'blue'))
                ax.XAxis.Color = 'b';
            else
                ax.XAxis.Color = [.6 0 0.4];
            end
            ah(j).handles = h;
            
            
            %     set(ah(j).ax, 'Children', h);
            
            
        end

    case '1stimRF_full_low_ker'
        
        SF = 60;
        
        
        for j = 1:length(btdstruct)
            
            
            if( strncmpi(btdstruct(j).var(StimNum).gqname, 'led2', 4))
                stim = 'blue';
            elseif( strncmpi(btdstruct(j).var(StimNum).gqname, 'led1', 4))
                stim = 'red';
            end
            
            ah(j).pos = pos{j}; %#ok<*AGROW>
            ah(j).ax = axes('Position', pos{j});
            ah(j).name = names{j};
            
            
            if(isfield(btdstruct(j), 'RatesU'))
                Rates = btdstruct(j).RatesU;
                stim = 'purple';
            elseif(isfield(btdstruct(j), 'Rates'))
                Rates = btdstruct(j).Rates;
            else
                Rates = btdstruct(j);
            end
            
            if(nargin>5)
                k = StimNum;
            else
                k=1;
            end
            
            lx_low_full = Rates.rfLowVar(k).lx;
            lx_low_low = Rates.rfLowVar_LK(k).lx;
            rate_low_full = SF*Rates.rfLowVar(k).rate;
            rate_eb_low_full = SF*Rates.rfLowVar(k).rate_eb;
            rate_low_low = SF*Rates.rfLowVar_LK(k).rate;
            rate_eb_low_low = SF*Rates.rfLowVar_LK(k).rate_eb;
            
           
            if(pdegree==1)                
                fitRate_low_full = SF*exp(polyval(Rates.rfLowVar(k).rateFitLin, lx_low_full));            
                fitRate_low_low = SF*exp(polyval(Rates.rfLowVar_LK(k).rateFitLin, lx_low_low));
            elseif(pdegree==2)
                fitRate_low_full = SF*exp(polyval(Rates.rfLowVar(k).rateFitQuad, lx_low_full));
                fitRate_low_low = SF*exp(polyval(Rates.rfLowVar_LK(k).rateFitQuad, lx_low_low));
            end
            
            
            
            h(1) = errorbar (lx_low_full, rate_low_full, rate_eb_low_full, 'color', colors{j}{1}, 'marker', 'o', 'linestyle', 'none'); %'MarkerSize', 5, 'CapSize', 5
            hold all;
            h(2) = errorbar (lx_low_low, rate_low_low, rate_eb_low_low, 'color', colors{j}{2}, 'marker', 's', 'MarkerFaceColor', colors{j}{2}, 'linestyle', 'none');
            
            if(pdegree>0)
                h(3) = plot (lx_low_full, fitRate_low_full, 'color', colors{j}{1}, 'linestyle', ':', 'LineWidth', 2);
                h(4) = plot (lx_low_full, fitRate_low_low, 'color', colors{j}{2}, 'linestyle', '-', 'LineWidth', 2);
            end
            
            title(names{j});
            
            ax = gca;
            ax.Box = 'off';
            %             axis square
            
            ax.YLim = [0 15];
            %             ax.YTick = [0 5 10];
            ax.YLabel.String = 'turn rate (min^{-1})';
            ax.YLabel.FontSize = 10;
            ax.XLabel.String = 'Filtered Stimulus Value';
            ax.XLabel.FontSize = 10;
            
            if(strcmpi(stim, 'red'))
                ax.XAxis.Color = 'r';
            elseif(strcmpi(stim, 'blue'))
                ax.XAxis.Color = 'b';
            else
                ax.XAxis.Color = [.6 0 0.4];
            end
            ah(j).handles = h;
            
            
            %     set(ah(j).ax, 'Children', h);
            
            
        end

    case '1stimRF_full_high_ker'
        
        SF = 60;
        
        
        for j = 1:length(btdstruct)
            
            
            if( strncmpi(btdstruct(j).var(StimNum).gqname, 'led2', 4))
                stim = 'blue';
            elseif( strncmpi(btdstruct(j).var(StimNum).gqname, 'led1', 4))
                stim = 'red';
            end
            
            ah(j).pos = pos{j}; %#ok<*AGROW>
            ah(j).ax = axes('Position', pos{j});
            ah(j).name = names{j};
            
            
            if(isfield(btdstruct(j), 'RatesU'))
                Rates = btdstruct(j).RatesU;
                stim = 'purple';
            elseif(isfield(btdstruct(j), 'Rates'))
                Rates = btdstruct(j).Rates;
            else
                Rates = btdstruct(j);
            end
            
            if(nargin>5)
                k = StimNum;
            else
                k=1;
            end
            
            
            lx_high_full = Rates.rfHighVar(k).lx;
            lx_high_high = Rates.rfHighVar_HK(k).lx;
            rate_high_full = SF*Rates.rfHighVar(k).rate;
            rate_eb_high_full = SF*Rates.rfHighVar(k).rate_eb;
            rate_high_high = SF*Rates.rfHighVar_HK(k).rate;
            rate_eb_high_high = SF*Rates.rfHighVar_HK(k).rate_eb;
            
            if(pdegree==1)
                
                fitRate_high_full = SF*exp(polyval(Rates.rfHighVar(k).rateFitLin, lx_high_full));
                fitRate_high_high = SF*exp(polyval(Rates.rfHighVar_HK(k).rateFitLin, lx_high_high));
            elseif(pdegree==2)
                fitRate_high_full = SF*exp(polyval(Rates.rfHighVar(k).rateFitQuad, lx_high_full));
                fitRate_high_high = SF*exp(polyval(Rates.rfHighVar_HK(k).rateFitQuad, lx_high_high));
            end
            
            
            
            h(1) = errorbar (lx_high_full, rate_high_full, rate_eb_high_full, 'color', colors{j}{1}, 'marker', 'o', 'linestyle', 'none'); %'MarkerSize', 5, 'CapSize', 5
            hold all;
            h(2) = errorbar (lx_high_high, rate_high_high, rate_eb_high_high, 'color', colors{j}{2}, 'marker', 's', 'MarkerFaceColor', colors{j}{2}, 'linestyle', 'none');
            
            if(pdegree>0)
                h(3) = plot (lx_high_full, fitRate_high_full, 'color', colors{j}{1}, 'linestyle', ':', 'LineWidth', 2);
                h(4) = plot (lx_high_full, fitRate_high_high, 'color', colors{j}{2}, 'linestyle', '-', 'LineWidth', 2);
            end
            
            title(names{j});
            
            ax = gca;
            ax.Box = 'off';
            %             axis square
            
            ax.YLim = [0 15];
            %             ax.YTick = [0 5 10];
            ax.YLabel.String = 'turn rate (min^{-1})';
            ax.YLabel.FontSize = 10;
            ax.XLabel.String = 'Filtered Stimulus Value';
            ax.XLabel.FontSize = 10;
            
            if(strcmpi(stim, 'red'))
                ax.XAxis.Color = 'r';
            elseif(strcmpi(stim, 'blue'))
                ax.XAxis.Color = 'b';
            else
                ax.XAxis.Color = [.6 0 0.4];
            end
            ah(j).handles = h;
            
            
            %     set(ah(j).ax, 'Children', h);
            
            
        end        

    case '1stimRF_both_ker'
        
        SF = 60;
        
        
        for j = 1:length(btdstruct)
            
            
            if( strncmpi(btdstruct(j).var(StimNum).gqname, 'led2', 4))
                stim = 'blue';
            elseif( strncmpi(btdstruct(j).var(StimNum).gqname, 'led1', 4))
                stim = 'red';
            end
            
            ah(j).pos = pos{j}; %#ok<*AGROW>
            ah(j).ax = axes('Position', pos{j});
            ah(j).name = names{j};
            
            
            if(isfield(btdstruct(j), 'RatesU'))
                Rates = btdstruct(j).RatesU;
                stim = 'purple';
            elseif(isfield(btdstruct(j), 'Rates'))
                Rates = btdstruct(j).Rates;
            else
                Rates = btdstruct(j);
            end
            
            if(nargin>5)
                k = StimNum;
            else
                k=1;
            end
            
            
            lx_high_full = Rates.rfHighVar(k).lx;
            lx_high_high = Rates.rfHighVar_HK(k).lx;
            rate_high_full = SF*Rates.rfHighVar(k).rate;
            rate_eb_high_full = SF*Rates.rfHighVar(k).rate_eb;
            rate_high_high = SF*Rates.rfHighVar_HK(k).rate;
            rate_eb_high_high = SF*Rates.rfHighVar_HK(k).rate_eb;
            
            lx_low_full = Rates.rfLowVar(k).lx;
            lx_low_low = Rates.rfLowVar_LK(k).lx;
            rate_low_full = SF*Rates.rfLowVar(k).rate;
            rate_eb_low_full = SF*Rates.rfLowVar(k).rate_eb;
            rate_low_low = SF*Rates.rfLowVar_LK(k).rate;
            rate_eb_low_low = SF*Rates.rfLowVar_LK(k).rate_eb;
            
            if(pdegree==1)
                
                fitRate_high_full = SF*exp(polyval(Rates.rfHighVar(k).rateFitLin, lx_high_full));
                fitRate_high_high = SF*exp(polyval(Rates.rfHighVar_HK(k).rateFitLin, lx_high_high));
                fitRate_low_full = SF*exp(polyval(Rates.rfLowVar(k).rateFitLin, lx_low_full));            
                fitRate_low_low = SF*exp(polyval(Rates.rfLowVar_LK(k).rateFitLin, lx_low_low));
            elseif(pdegree==2)
                fitRate_high_full = SF*exp(polyval(Rates.rfHighVar(k).rateFitQuad, lx_high_full));
                fitRate_high_high = SF*exp(polyval(Rates.rfHighVar_HK(k).rateFitQuad, lx_high_high));
                fitRate_low_full = SF*exp(polyval(Rates.rfLowVar(k).rateFitQuad, lx_low_full));
                fitRate_low_low = SF*exp(polyval(Rates.rfLowVar_LK(k).rateFitQuad, lx_low_low));
            end
            
            
            
            h(1) = errorbar (lx_high_full, rate_high_full, rate_eb_high_full, 'color', colors{1}, 'marker', 'o', 'linestyle', 'none'); %'MarkerSize', 5, 'CapSize', 5
            hold all;
            h(2) = errorbar (lx_high_high, rate_high_high, rate_eb_high_high, 'color', colors{2}, 'marker', 's', 'MarkerFaceColor', colors{2}, 'linestyle', 'none');
            h(3) = errorbar (lx_low_full, rate_low_full, rate_eb_low_full, 'color', colors{3}, 'marker', 'o', 'linestyle', 'none'); %'MarkerSize', 5, 'CapSize', 5
            h(4) = errorbar (lx_low_low, rate_low_low, rate_eb_low_low, 'color', colors{4}, 'marker', 's', 'MarkerFaceColor', colors{4}, 'linestyle', 'none');
            if(pdegree>0)
                h(5) = plot (lx_high_full, fitRate_high_full, 'color', colors{1}, 'linestyle', ':', 'LineWidth', 2);
                h(6) = plot (lx_high_full, fitRate_high_high, 'color', colors{2}, 'linestyle', '-', 'LineWidth', 2);
                h(7) = plot (lx_low_full, fitRate_low_full, 'color', colors{3}, 'linestyle', ':', 'LineWidth', 2);
                h(8) = plot (lx_low_full, fitRate_low_low, 'color', colors{4}, 'linestyle', '-', 'LineWidth', 2);
            end
            
            title(names{j});
            
            ax = gca;
            ax.Box = 'off';
            %             axis square
            
            ax.YLim = [0 15];
            %             ax.YTick = [0 5 10];
            ax.YLabel.String = 'turn rate (min^{-1})';
            ax.YLabel.FontSize = 10;
            ax.XLabel.String = 'Filtered Stimulus Value';
            ax.XLabel.FontSize = 10;
            
            if(strcmpi(stim, 'red'))
                ax.XAxis.Color = 'r';
            elseif(strcmpi(stim, 'blue'))
                ax.XAxis.Color = 'b';
            else
                ax.XAxis.Color = [.6 0 0.4];
            end
            ah(j).handles = h;
            
            
            %     set(ah(j).ax, 'Children', h);
            
            
        end
        
    case '1stimRF_levels'
        
        SF = 60;
        
        
        for j = 1:length(btdstruct)
            
            
            if( strncmpi(btdstruct(j).var(StimNum).gqname, 'led2', 4))
                stim = 'blue';
            elseif( strncmpi(btdstruct(j).var(StimNum).gqname, 'led1', 4))
                stim = 'red';
            end
            
            ah(j).pos = pos{j}; %#ok<*AGROW>
            ah(j).ax = axes('Position', pos{j});
            ah(j).name = names{j};
            
            
            if(isfield(btdstruct(j), 'RatesU'))
                Rates = btdstruct(j).RatesU;
                stim = 'purple';
            elseif(isfield(btdstruct(j), 'Rates'))
                Rates = btdstruct(j).Rates;
            else
                Rates = btdstruct(j);
            end
            
            if(nargin>5)
                k = StimNum;
            else
                k=1;
            end
            
            lx_low = Rates.rfLowVar_lvl(k).lx;
            lx_high = Rates.rfHighVar_lvl(k).lx;
            rate_low = SF*Rates.rfLowVar_lvl(k).rate;
            rate_eb_low = SF*Rates.rfLowVar_lvl(k).rate_eb;
            rate_high = SF*Rates.rfHighVar_lvl(k).rate;
            rate_eb_high = SF*Rates.rfHighVar_lvl(k).rate_eb;
            
            if(pdegree==1)
                
                fitRate_low = SF*exp(polyval(Rates.rfLowVar_lvl(k).rateFitLin, lx_low));
                fitRate_high = SF*exp(polyval(Rates.rfHighVar_lvl(k).rateFitLin, lx_high));
            elseif(pdegree==2)
                fitRate_low = SF*exp(polyval(Rates.rfLowVar_lvl(k).rateFitQuad, lx_low));
                fitRate_high = SF*exp(polyval(Rates.rfHighVar_lvl(k).rateFitQuad, lx_high));
            end
            
            h(1) = errorbar (lx_low, rate_low, rate_eb_low, 'color', colors{j}{1}, 'marker', 'o', 'linestyle', 'none'); %'MarkerSize', 5, 'CapSize', 5
            hold all;
            h(2) = errorbar (lx_high, rate_high, rate_eb_high, 'color', colors{j}{2}, 'marker', 's', 'MarkerFaceColor', colors{j}{2}, 'linestyle', 'none');
            
            if(pdegree>0)
                h(3) = plot (lx_low, fitRate_low, 'color', colors{j}{1}, 'linestyle', ':', 'LineWidth', 2);
                h(4) = plot (lx_high, fitRate_high, 'color', colors{j}{2}, 'linestyle', '-', 'LineWidth', 2);
            end
            
            title(names{j});
            
            ax = gca;
            ax.Box = 'off';
            %             axis square
            
            ax.YLim = [0 15];
            %             ax.YTick = [0 5 10];
            ax.YLabel.String = 'turn rate (min^{-1})';
            ax.YLabel.FontSize = 10;
            ax.XLabel.String = 'Filtered Stimulus Value';
            ax.XLabel.FontSize = 10;
            
            if(strcmpi(stim, 'red'))
                ax.XAxis.Color = 'r';
            elseif(strcmpi(stim, 'blue'))
                ax.XAxis.Color = 'b';
            else
                ax.XAxis.Color = [.6 0 0.4];
            end
            ah(j).handles = h;
            
            
            %     set(ah(j).ax, 'Children', h);
            
            
        end        

    case '1stimRF_lev_120kernel'
        
        SF = 60;
        
        
        for j = 1:length(btdstruct)
            
            
            if( strncmpi(btdstruct(j).var(StimNum).gqname, 'led2', 4))
                stim = 'blue';
            elseif( strncmpi(btdstruct(j).var(StimNum).gqname, 'led1', 4))
                stim = 'red';
            end
            
            ah(j).pos = pos{j}; %#ok<*AGROW>
            ah(j).ax = axes('Position', pos{j});
            ah(j).name = names{j};
            
            
            if(isfield(btdstruct(j), 'RatesU'))
                Rates = btdstruct(j).RatesU;
                stim = 'purple';
            elseif(isfield(btdstruct(j), 'Rates'))
                Rates = btdstruct(j).Rates;
            else
                Rates = btdstruct(j);
            end
            
            if(nargin>5)
                k = StimNum;
            else
                k=1;
            end
            
            lx_low = Rates.rfLowVar_ber120(k).lx;
            lx_high = Rates.rfHighVar_ber120(k).lx;
            rate_low = SF*Rates.rfLowVar_ber120(k).rate;
            rate_eb_low = SF*Rates.rfLowVar_ber120(k).rate_eb;
            rate_high = SF*Rates.rfHighVar_ber120(k).rate;
            rate_eb_high = SF*Rates.rfHighVar_ber120(k).rate_eb;
            
            if(pdegree==1)
                
                fitRate_low = SF*exp(polyval(Rates.rfLowVar(k).rateFitLin, lx_low));
                fitRate_high = SF*exp(polyval(Rates.rfHighVar(k).rateFitLin, lx_high));
            elseif(pdegree==2)
                fitRate_low = SF*exp(polyval(Rates.rfLowVar(k).rateFitQuad, lx_low));
                fitRate_high = SF*exp(polyval(Rates.rfHighVar(k).rateFitQuad, lx_high));
            end
            
            h(1) = errorbar (lx_low, rate_low, rate_eb_low, 'color', colors{j}{1}, 'marker', 'o', 'linestyle', 'none'); %'MarkerSize', 5, 'CapSize', 5
            hold all;
            h(2) = errorbar (lx_high, rate_high, rate_eb_high, 'color', colors{j}{2}, 'marker', 's', 'MarkerFaceColor', colors{j}{2}, 'linestyle', 'none');
            
            if(pdegree>0)
                h(3) = plot (lx_low, fitRate_low, 'color', colors{j}{1}, 'linestyle', ':', 'LineWidth', 2);
                h(4) = plot (lx_high, fitRate_high, 'color', colors{j}{2}, 'linestyle', '-', 'LineWidth', 2);
            end
            
            title(names{j});
            
            ax = gca;
            ax.Box = 'off';
            %             axis square
            
            ax.YLim = [0 15];
            %             ax.YTick = [0 5 10];
            ax.YLabel.String = 'turn rate (min^{-1})';
            ax.YLabel.FontSize = 10;
            ax.XLabel.String = 'Filtered Stimulus Value';
            ax.XLabel.FontSize = 10;
            
            if(strcmpi(stim, 'red'))
                ax.XAxis.Color = 'r';
            elseif(strcmpi(stim, 'blue'))
                ax.XAxis.Color = 'b';
            else
                ax.XAxis.Color = [.6 0 0.4];
            end
            ah(j).handles = h;
            
            
            %     set(ah(j).ax, 'Children', h);
            
            
        end        
        
    case '1stimScaledRF'
        
        for j=1:length(btdstruct)
            
            SF = 60*deltaT/0.05;
            
            ah(j).pos = pos{j}; %#ok<*AGROW>
            ah(j).ax = axes('Position', ah(j).pos);
            ah(j).name = names{j};
            
            if(isfield(btdstruct(j), 'ScaledRatesU'))
                ScaledRates = btdstruct(j).ScaledRatesU;
            elseif(isfield(btdstruct(j), 'ScaledRates'))
                ScaledRates = btdstruct(j).ScaledRates;
            else
                ScaledRates = btdstruct(j);
            end
            rfH = ScaledRates.scaledRFHigh(1);
            rfL = ScaledRates.scaledRFLow(1);
            rfU = ScaledRates.scaledRFUp(1);
            rfD = ScaledRates.scaledRFDown(1);
            
            if(pdegree == 1)
                fitH = rfH.rateFitLin;
                fitL = rfL.rateFitLin;
                fitU = rfU.rateFitLin;
                fitD = rfD.rateFitLin;
            else
                fitH = rfH.rateFitQuad;
                fitL = rfL.rateFitQuad;
                fitU = rfU.rateFitQuad;
                fitD = rfD.rateFitQuad;
            end
            
            h(1) = errorbar(rfH.lx, SF*rfH.rate, SF*rfH.rate_eb, 'color', colors{j}{1}, 'marker', 'o', 'linestyle', 'none');
            hold all
            h(2) = errorbar(rfL.lx, SF*rfL.rate, SF*rfL.rate_eb, 'color', colors{j}{2}, 'marker', 's', 'linestyle', 'none');
            h(3) = errorbar(rfU.lx, SF*rfU.rate, SF*rfU.rate_eb, 'color', 'm', 'marker', '^', 'linestyle', 'none');
            h(4) = errorbar(rfD.lx, SF*rfD.rate, SF*rfD.rate_eb, 'color', 'c', 'marker', 'v', 'linestyle', 'none');
            %
            h(5) = plot(rfH.lx, SF*exp(polyval(fitH, rfH.lx)), 'color', colors{j}{1});
            h(6) = plot( rfL.lx, SF*exp(polyval(fitL, rfL.lx)), 'color', colors{j}{2});
            h(7) = plot(rfU.lx, SF*exp(polyval(fitU, rfU.lx)), 'color', 'm');
            h(8) = plot( rfD.lx, SF*exp(polyval(fitD, rfD.lx)), 'color', 'c');
            
            ax = gca;
            %
            %             if(log)
            %                 ax.YScale = 'log';
            %                 ax.YLim = [0.4 20];
            %                 ax.YTick = [1 2 5 10 20];
            %                 ax.YMinorTick = 'off';
            %                 ax.YTickLabel = {'1', '2', '5', '10', '20', '30'};
            %             end
            ax.YLabel.String = 'turn rate (min^{-1})';
            ax.YLabel.FontSize = 10;
            ax.XLabel.String = 'Scaled Stimulus Value - \alpha\cdotx';
            ax.XLabel.FontSize = 10;
            ax.Box = 'off';
            
            
            title(names{j});
            ah(j).handles = h;
        end
        
    case '2stimPPF'
        
        for j=1:length(btdstruct)
            
            tshift = 3;
            
            period = btdstruct(j).var.period;
            
            if(strcmpi(estimator, 'PPF'))
                tx = btdstruct(j).alphavseti.tx_ton;
                dt = round(period)/length(tx);
                alpha = btdstruct(j).alphavseti.alpha_ton(1,:);
                beta = btdstruct(j).alphavseti.alpha_ton(2,:);
                valpha = sqrt(squeeze(btdstruct(j).alphavseti.valpha_ton(1,1,:)));
                vbeta = sqrt(squeeze(btdstruct(j).alphavseti.valpha_ton(2,2,:)));
            elseif(strcmpi(estimator, 'Bayes'))
                if(strncmpi(multiType, 'add', 3))
                    BayesStruct = btdstruct(j).BayesAlphavsEti_OL.add;
                elseif(strncmpi(multiType, 'mult', 3))
                    BayesStruct = btdstruct(j).BayesAlphavsEti_OL.multiply;
                else
                    BayesStruct = btdstruct(j).BayesAlphavsEti_OL;
                end
                tx = BayesStruct.tx_ton;
                dt = round(period)/length(tx);
                alpha = BayesStruct.alpha_ton(1,:);
                beta = BayesStruct.alpha_ton(2,:);
                valpha = sqrt(squeeze(BayesStruct.valpha_ton(1,1,:)))';
                vbeta = sqrt(squeeze(BayesStruct.valpha_ton(2,2,:)))';
                
            end

            tx_shift = -fliplr(tx);
            %             tEst = [-tshift:dt:tx_ton(1), tx_ton];
            alpha_shift = circshift(alpha, round(length(tx)/2));
            beta_shift = circshift(beta, round(length(tx)/2));
            valpha_shift = circshift(valpha, round(length(tx)/2));
            vbeta_shift = circshift(vbeta, round(length(tx)/2));
            
            tEst = [tx_shift(round(end-tshift/dt):end), tx];
            aEst = [alpha_shift(round(end-tshift/dt):end), alpha_shift];
            vaEst = [valpha_shift(round(end-tshift/dt):end), valpha_shift];
            bEst = [beta_shift(round(end-tshift/dt):end), beta_shift];
            vbEst = [vbeta_shift(round(end-tshift/dt):end), vbeta_shift];

            
            ah(j).pos = pos{j}; %#ok<*AGROW>
            ah(j).ax = axes('Position', ah(j).pos);
            ah(j).name = names{j};
            
            h(1,:) = shadedErrorPlot({tEst tEst}, {aEst bEst}, {vaEst, vbEst}, {vaEst, vbEst}, {colors{1}, colors{2}});
            
            
            pp = h(1,:).Parent;
            pp.XTick = [0 20 pp.XLim(2)];
            pp.XTickLabel = {'0', '20', num2str(round(period))};
            
            pp.Children(1).LineWidth = 1.5;
            pp.Children(2).LineWidth = 1.5;
            
            ax = gca;
            ax.YLim = [0.6 1.7];
            hold on
            plot([0 0], [ax.YLim(1) ax.YLim(2)], 'k--', 'linewidth', 1.5);
            plot([period/2 period/2], [ax.YLim(1) ax.YLim(2)], 'k--', 'linewidth', 1.5);
            
            ax.YLabel.String = '\alpha';
            ax.YLabel.FontSize = 10;
            ax.XLabel.String = 'time in cycle (s)';
            ax.XLabel.FontSize = 10;
            
            
            ax.DataAspectRatio = [30 1 1];
            
            ax.Box = 'off';
            title(names{j});
            ah(j).handles = h;
            
            
        end
        
        
        
    case 'TTA'
        
        
        lowp = @(x) lowpass1D(x, 2.5);
        
        for j=1:length(btdstruct)
            
                    
            if( strncmpi(btdstruct(j).var.gqname, 'led2', 4))
                stim = 'blue';
            elseif( strncmpi(btdstruct(j).var.gqname, 'led1', 4))
                stim = 'red';
            end
            
            ah(j).pos = pos{j}; %#ok<*AGROW>
            ah(j).ax = axes('Position', pos{j});
            ah(j).name = names{j};
            
            ttaxis = btdstruct(j).Kernels.High.taxis;
            ttaL = btdstruct(j).Kernels.Low.tta;
            ttaH = btdstruct(j).Kernels.High.tta;
            SFLmax = max(abs(lowp(ttaL)));
            SFHmax = max(abs(lowp(ttaH)));
            SFLmin = min(lowp(ttaL));
            SFHmin = min(lowp(ttaH));
            ttaL = ttaL/SFLmax;
            ttaH = ttaH/SFHmax;

            
            
            h(1) = plot (ttaxis, lowp(ttaL), 'color', colors{j}{1}, 'LineWidth', 1.5);
            hold all
            h(2) = plot (ttaxis, lowp(ttaH), 'color', colors{j}{2}, 'LineWidth', 1.5);
            
            
            title(names{j});
            
            %             axis square
            ax = gca;
            ax.Box = 'off';
            ax.XLim = [min(ttaxis) max(ttaxis)];
            ax.XLim = [-6 1];
            ax.YLim = [min(ttaH) max(ttaH)];
            ax.YTick = 0;
            ax.YTickLabel = '0';
            
            if(strcmpi(stim, 'red'))
                ax.YAxis.Color = 'r';
                ax.YLabel.String = {'Avg. Red'; 'Intensity Derivative (Scaled)'};
            elseif(strcmpi(stim, 'blue'))
                ax.YAxis.Color = 'b';
                ax.YLabel.String = {'Avg. Blue'; 'Intensity Derivative (Scaled)'};
                if(strcmpi(btdstruct(j).var.gqname, 'led2Val'))
                    ax.YAxis.Color = 'b';
                    ax.YLabel.String = {'Avg. Blue'; 'Intensity (Scaled)'};
                end
            end
            
            ax.YLabel.FontSize = 8;
            ax.XLabel.String = {'time (s)'; 'rel. to turn start'};
            ax.XLabel.FontSize = 10;
            h(3) = plot([0 0], [min(ax.YLim) max(ax.YLim)], 'k:', 'LineWidth', 1.5);
            h(4) = plot([min(ax.XLim) max(ax.XLim)], [0 0], 'k--', 'LineWidth', 0.1);
            
            ah(j).handles = h;
            
        end
        
    case 'TTA_levels'
        
        
        lowp = @(x) lowpass1D(x, 2.5);
        
        for j=1:length(btdstruct)
            
                    
            if( strncmpi(btdstruct(j).var.gqname, 'led2', 4))
                stim = 'blue';
            elseif( strncmpi(btdstruct(j).var.gqname, 'led1', 4))
                stim = 'red';
            end
            
            ah(j).pos = pos{j}; %#ok<*AGROW>
            ah(j).ax = axes('Position', pos{j});
            ah(j).name = names{j};
            
            ttaxis = btdstruct(j).Kernels.High_lvl.taxis;
            ttaL = btdstruct(j).Kernels.Low_lvl.tta;
            ttaH = btdstruct(j).Kernels.High_lvl.tta;
            SFLmax = max(lowp(ttaL));
            SFHmax = max(lowp(ttaH));
            SFLmin = min(lowp(ttaL));
            SFHmin = min(lowp(ttaH));
            ttaL = ttaL/(SFLmax-SFLmin);
            ttaH = ttaH/(SFHmax-SFHmin);
            
            h(1) = plot (ttaxis, lowp(ttaL), 'color', colors{j}{1}, 'LineWidth', 1.5);
            hold all
            h(2) = plot (ttaxis, lowp(ttaH), 'color', colors{j}{2}, 'LineWidth', 1.5);
            
            
            title(names{j});
            
            %             axis square
            ax = gca;
            ax.Box = 'off';
            ax.XLim = [min(ttaxis) max(ttaxis)];
            ax.XLim = [-6 1];
            ax.YLim = [min(ttaH) max(ttaH)];
            ax.YTick = 0;
            ax.YTickLabel = '0';
            
            if(strcmpi(stim, 'red'))
                ax.YAxis.Color = 'r';
                ax.YLabel.String = {'Avg. Red'; 'Intensity (Scaled)'};
            elseif(strcmpi(stim, 'blue'))
                ax.YAxis.Color = 'b';
                ax.YLabel.String = {'Avg. Blue'; 'Intensity (Scaled)'};
            end
            
            ax.YLabel.FontSize = 8;
            ax.XLabel.String = {'time (s)'; 'rel. to turn start'};
            ax.XLabel.FontSize = 10;
            h(3) = plot([0 0], [min(ax.YLim) max(ax.YLim)], 'k:', 'LineWidth', 1.5);
            h(4) = plot([min(ax.XLim) max(ax.XLim)], [0 0], 'k:', 'LineWidth', 1.5);
            
            ah(j).handles = h;
            
        end
        
    case 'TTA_bootstrap'
        
        
        lowp = @(x) lowpass1D(x, 2.5);
        
            
            
            ah.pos = pos{1}; %#ok<*AGROW>
            ah.ax = axes('Position', pos{1});
            ah.name = names{1};
            
            ttaxis = btdstruct.High.taxis;
            ttaL = btdstruct.Low.bootstrapMean;
            ttaH = btdstruct.High.bootstrapMean;
            ttaLs = btdstruct.Low.bootstrapStd;
            ttaHs = btdstruct.High.bootstrapStd;
            SFLmax = max(abs(lowp(ttaL)));
            SFHmax = max(abs(lowp(ttaH)));
            SFLmin = min(lowp(ttaL));
            SFHmin = min(lowp(ttaH));
            ttaL = ttaL/(SFLmax);
            ttaH = ttaH/(SFHmax);
            ttaLs = ttaLs/(SFLmax);
            ttaHs = ttaHs/(SFHmax);
            
            temp = shadedErrorPlot({ttaxis ttaxis}, {lowp(ttaL) lowp(ttaH)}, {lowp(ttaLs) lowp(ttaHs)});
            h(1) = temp(1);
            h(2) = temp(2);
            h(3) = temp(3);
            
            title(names{1});
            pp = h(1).Parent;
%                pp.Children(1).Color = 1/255*[74 137 165];
%                pp.Children(2).Color = [.8 0.2 0.6];
%                pp.Children(1).LineWidth = 1.5;
%                pp.Children(2).LineWidth = 1.2;
%                pp.Children(3).FaceColor = 1/255*[57 70 156];
%                pp.Children(4).FaceColor = [.4 0 0.6];
%                pp.Children(3).FaceAlpha = 0.5;
%                pp.Children(4).FaceAlpha = .3;
                pp.Children(1).Color = colors{1}{1};
                pp.Children(2).Color = colors{1}{2};
                pp.Children(1).LineWidth = 1.5;
                pp.Children(2).LineWidth = 1.2;
                pp.Children(3).FaceColor = colors{1}{1};
                pp.Children(4).FaceColor = colors{1}{2};
                pp.Children(3).FaceAlpha = 0.5;
                pp.Children(4).FaceAlpha = 0.3;
            %             axis square
            ax = gca;
            ax.Box = 'off';
            ax.XLim = [min(ttaxis) max(ttaxis)];
            ax.XLim = [-6 1];
            ax.YLim = [min(ttaH) max(ttaH)];
            ax.YTick = 0;
            ax.YTickLabel = '0';
            
            if(strncmpi(names, 'Berlin',6))
                ax.YAxis.Color = 'b';
                ax.YLabel.String = {'Avg. Blue'; 'Intensity Derivative (Scaled)'};
            else
                ax.YAxis.Color = 'r';
                ax.YLabel.String = {'Avg. Red'; 'Intensity Derivative (Scaled)'};
            end                      
                        
            ax.YLabel.FontSize = 8;
            ax.XLabel.String = {'time (s)'; 'rel. to turn start'};
            ax.XLabel.FontSize = 10;
            hold on;
            h(4) = plot([0 0], [min(ax.YLim) max(ax.YLim)], 'k:', 'LineWidth', 1.5);
            h(5) = plot([min(ax.XLim) max(ax.XLim)], [0 0], 'k--', 'LineWidth', 0.1);
            
            ah.handles = h;
            hold off;
       
        
    case 'TTA-AccRej'
        
        for j=1:length(btdstruct)
            ah(j).pos = pos{j}; %#ok<*AGROW>
            ah(j).ax = axes('Position', pos{j});
            ah(j).name = names{j};
            
            taxis = btdstruct(j).Kernels.Low.taxis;
            ttaL = btdstruct(j).Kernels.Low.tta;
            ttaH = btdstruct(j).Kernels.High.tta;
            ttaD = btdstruct(j).Kernels.Down.tta;
            ttaU = btdstruct(j).Kernels.Up.tta;
            
            
            AccL = btdstruct(j).Kernels.Low.tta_acchs; AccL = AccL/max(abs(AccL));
            RejL = btdstruct(j).Kernels.Low.tta_rejhs; RejL = RejL/max(abs(RejL));
            AccH = btdstruct(j).Kernels.High.tta_acchs; AccH = AccH/max(abs(AccH));
            RejH = btdstruct(j).Kernels.High.tta_rejhs; RejH = RejH/max(abs(RejH));
            AccU = btdstruct(j).Kernels.Up.tta_acchs; AccU = AccU/max(abs(AccU));
            RejU = btdstruct(j).Kernels.Up.tta_rejhs; RejU = RejU/max(abs(RejU));
            AccD = btdstruct(j).Kernels.Down.tta_acchs; AccD = AccD/max(abs(AccD));
            RejD = btdstruct(j).Kernels.Down.tta_rejhs; RejD = RejD/max(abs(RejD));
            
            %             SFL = max( max(abs(AccL)), max(abs(RejL)) );
            %             SFH = max( max(abs(AccH)), max(abs(RejH)) );
            %             SFU = max( max(abs(AccU)), max(abs(RejU)) );
            %             SFD = max( max(abs(AccD)), max(abs(RejD)) );
            %             AccL = AccL/SFL;
            %             RejL = RejL/SFL;
            %             AccH = AccH/SFH;
            %             RejH = RejH/SFH;
            %             AccU = AccU/SFU;
            %             RejU = RejU/SFU;
            %             AccD = AccD/SFD;
            %             RejD = RejD/SFD;
            
            h(1) = plot (taxis, lowpass1D(AccL, 2.5), ':', 'color', colors{j}{1}, 'LineWidth', 1.5);
            hold all
            h(2) = plot (taxis, lowpass1D(RejL, 2.5), '--', 'color', colors{j}{1}, 'LineWidth', 1.5);
            h(3) = plot (taxis, lowpass1D(AccH, 2.5), ':', 'color', colors{j}{2}, 'LineWidth', 1.5);
            h(4) = plot (taxis, lowpass1D(RejH, 2.5), '--', 'color', colors{j}{2}, 'LineWidth', 1.5);
            
            ax = gca;
            ax.Box = 'off';
            ax.XLim = [-7 6];
            ax.YLim = [min([min(RejH), min(AccH), min(AccL), min(RejL)]) max([max(RejH), max(AccH), max(AccL), max(RejL)])];
            ax.YTick = 0;
            ax.YTickLabel = '0';
            ax.YAxis.Color = 'r';
            
            ax.YLabel.String = {'Avg. Red'; 'Intensity Derivative (Scaled)'};
            ax.YLabel.FontSize = 8;
            ax.XLabel.String = {'time (s)'; 'rel. to turn start'};
            ax.XLabel.FontSize = 10;
            h(3) = plot([0 0], [min(ax.YLim) max(ax.YLim)], 'k:', 'LineWidth', 1.5);
            h(4) = plot([min(ax.XLim) max(ax.XLim)], [0 0], 'k--', 'LineWidth', 0.1);
            
            ah(j).handles = h;
            
            title(names{j});
            
        end
        
    case 'VarRamp'
        
        for j = 1:length(btdstruct)
            
            clear h h2 h3
            nstim = size(btdstruct(j).var, 2);
            
            SF = 60;
            
            k = j+length(btdstruct);
            kk = j+2*length(btdstruct);
            
            ah(j).pos = pos{j}; %#ok<*AGROW>
            ah(j).ax = axes('Position', pos{j});
            ah(j).name = names{j};
            
            
            
            if(strcmpi(estimator, 'PPF'))
                GR = btdstruct(j).gainrescaling_PPF;
                if(nstim==1)
                    ton = btdstruct(j).alphavseti.tx_ton;
                    alpha_s = btdstruct(j).alphavseti.alpha_s_ton;
                    valpha_s = sqrt(squeeze(btdstruct(j).alphavseti.valpha_s_ton));
                else
                    ton = btdstruct(j).alphavseti_U.tx_ton;
                    alpha_s = btdstruct(j).alphavseti_U.alpha_s_ton;
                    valpha_s = sqrt(squeeze(btdstruct(j).alphavseti_U.valpha_s_ton));
                end
            else
                GR = btdstruct(j).gainrescaling_Bayes;
                if(nstim==1)
                    ton = btdstruct(j).BayesAlphavsEti_OL.tx_ton;
                    alpha_s = btdstruct(j).BayesAlphavsEti_OL.alpha_ton;
                    valpha_s = sqrt(squeeze(btdstruct(j).BayesAlphavsEti_OL.valpha_ton)');
                else
                    ton = btdstruct(j).BayesAlphavsEti_U.tx_ton;
                    alpha_s = btdstruct(j).BayesAlphavsEti_U.alpha_ton;
                    valpha_s = sqrt(squeeze(btdstruct(j).BayesAlphavsEti_U.valpha_ton));
                end
            end
            
            sall = GR.son_all;
            sup = sall(1:round(end/2)-1);
            sdown = sall(round(end/2+1):end);
            alpha_s_up = alpha_s(1:round(end/2)-1);
            alpha_s_down = alpha_s(round(end/2)+1:end);
            fP = GR.params.all;
            fitfun = GR.fitfun;
            
            [xu, aup, ~, aup_eb] = meanyvsx(sup, alpha_s_up, linspace(min(sup), max(sup), 8));
            [xd, adown, ~, adown_eb] = meanyvsx(sdown, alpha_s_down, linspace(min(sdown), max(sdown), 8));
            [xsall, aall, ~, aall_eb] = meanyvsx(sall, alpha_s, linspace(min(sall), max(sall), 8));
            
            %{
            if(nstim==1)
                rates = btdstruct(j).Rates;
%                 scaled_rates = btdstruct(j).ScaledRates;
            else
                rates = btdstruct(j).RatesU;
%                 scaled_rates = btdstruct(j).ScaledRatesU;
            end
            xH = rates.rfHighVar.lx;
            rH = SF*rates.rfHighVar.rate;
            rH_eb = SF*rates.rfHighVar.rate_eb;
            xL = rates.rfLowVar.lx;
            rL = SF*rates.rfLowVar.rate;
            rL_eb = SF*rates.rfLowVar.rate_eb;
            
            sxH = scaled_rates.scaledRFHigh.lx;
            srH = SF*scaled_rates.scaledRFHigh.rate;
            srH_eb = SF*scaled_rates.scaledRFHigh.rate_eb;
            sxL = scaled_rates.scaledRFLow.lx;
            srL = SF*scaled_rates.scaledRFLow.rate;
            srL_eb = SF*scaled_rates.scaledRFLow.rate_eb;

            
            if(nstim==2)
                ratesBR = btdstruct(j).Rates;
                scaled_ratesBR = btdstruct(j).ScaledRates;
                
                xHo = ratesBR.rfHighVar(1).lx;
                rHo = SF*ratesBR.rfHighVar(1).rate;
                rHo_eb = SF*ratesBR.rfHighVar(1).rate_eb;
                xLo = ratesBR.rfLowVar(1).lx;
                rLo = SF*ratesBR.rfLowVar(1).rate;
                rLo_eb = SF*ratesBR.rfLowVar(1).rate_eb;
                
                sxHo = scaled_ratesBR.scaledRFHigh(1).lx;
                srHo = SF*scaled_ratesBR.scaledRFHigh(1).rate;
                srHo_eb = SF*scaled_ratesBR.scaledRFHigh(1).rate_eb;
                sxLo = scaled_ratesBR.scaledRFLow(1).lx;
                srLo = SF*scaled_ratesBR.scaledRFLow(1).rate;
                srLo_eb = SF*scaled_ratesBR.scaledRFLow(1).rate_eb;
                
                xHl = ratesBR.rfHighVar(2).lx;
                rHl = SF*ratesBR.rfHighVar(2).rate;
                rHl_eb = SF*ratesBR.rfHighVar(2).rate_eb;
                xLl = ratesBR.rfLowVar(2).lx;
                rLl = SF*ratesBR.rfLowVar(2).rate;
                rLl_eb = SF*ratesBR.rfLowVar(2).rate_eb;
                
                sxHl = scaled_ratesBR.scaledRFHigh(2).lx;
                srHl = SF*scaled_ratesBR.scaledRFHigh(2).rate;
                srHl_eb = SF*scaled_ratesBR.scaledRFHigh(2).rate_eb;
                sxLl = scaled_ratesBR.scaledRFLow(2).lx;
                srLl = SF*scaled_ratesBR.scaledRFLow(2).rate;
                srLl_eb = SF*scaled_ratesBR.scaledRFLow(2).rate_eb;
                
                
            end
            %}
            
            mlow = 'o';
            mhigh = '*';
            cup = [0 0.8 0.4];
            cdown = [0 0.4 0];
            reds = {1/255*[243 127 129], 1/255*[127 20 22]};
            blues = {1/255*[74 137 165], 1/255*[57 70 156]};
            purples = {[.8 0.2 0.6], [.4 0 0.6]};
            
            if(strcmpi(colors{j}, 'Red'))
                colorsHL = reds;
                calpha = {[1 0 0], [1 .8 .8]};
                
            elseif(strcmpi(colors{j}, 'Blue'))
                colorsHL = blues;
                calpha = {[0 0 1], [.8 .8 1]};
                
            else
                colorsHL = purples;
                calpha = {[.6 0 0.4], [.6 0.2 0.4]};
            end
            
            
            markers = {mlow, mhigh, mhigh, mlow};
            colors1 = {cup, cup, cdown, cdown};
            colorsBR = {reds{1} reds{2} blues{1} blues{2}};
            
            yyaxis left
            h(1,:) = plot(ton, sall, 'k');
            
            yyaxis right
            if(size(alpha_s, 1) == size(valpha_s, 2))
                valpha_s = valpha_s';
            end
            temp = shadedErrorPlot(ton, alpha_s,  valpha_s, valpha_s, calpha{1});
            h(2) = temp(1);
            h(3) = temp(2);
            
            pp = h(2).Parent;
            pp.Children(1).LineWidth = 1.5;
            
            title(names{j});
            ax = gca;
            ax.Box = 'off';
            ax.XLim = [min(ton) max(ton)];
            ax.YAxis(1).Limits = [0.5 3.5];
            ax.YAxis(1).Color = 'k';
            ax.YAxis(2).Limits = [0.5 1.8];
            ax.YAxis(2).TickValues = [0.5 1 1.5];
            ax.YAxis(2).Color = calpha{1};
            ax.YAxis(1).Label.String = '\sigma_{stim}';
            ax.YAxis(2).Label.String = '\alpha';
            ax.YAxis(2).Label.Rotation = 0;
            ax.YAxis(1).Label.FontSize = 10;
            ax.YAxis(2).Label.FontSize = 10;
            ax.XLabel.String = 'time in cycle (s)';
            ax.XLabel.FontSize = 10;
            
            %         set(ax, 'YScale', 'log');
            %         ax.YLim = [0 10];
            %         ax.YTick = [0 1 2 5 10];
            %         ax.YTickLabel = {'0', '1', '2', '5', '10'};
            
            ah(j).handles = h;
            
            ah(k).pos = pos{k}; %#ok<*AGROW>
            ah(k).ax = axes('Position', pos{k});
            ah(k).name = names{j};
            h2(1) = errorbar(xu, aup, aup_eb, 'color', cup, 'linestyle', 'none', 'marker', '^', 'MarkerSize', 6);
            hold on;
            h2(2) = errorbar(xd, adown, adown_eb, 'color', cdown, 'linestyle', 'none', 'marker', 'v', 'MarkerSize', 6);
            h2(3) = errorbar(xsall, aall, aall_eb, 'color','k', 'linestyle', 'none', 'marker', '.', 'MarkerSize', 12);
            h2(4) = plot(xsall, fitfun(fP, xsall), 'k-');
            
            title(names{j});
            ax = gca;
            ax.Box = 'off';
            ax.XLim = h(1).Parent.YAxis(1).Limits;
            ax.XLim = [0.8 3];
            ax.YLim = h(2).Parent.YLim;
            ax.YTick = h(2).Parent.YTick;
            ax.YAxis.Color = calpha{1};
            ax.XLabel.String = '\sigma_{stim}';
            ax.XLabel.FontSize = 10;
            
            ah(k).handles = h2;
            
            continue
            ah(kk).pos = pos{kk}; %#ok<*AGROW>
            ah(kk).ax = axes('Position', pos{kk});
            ah(kk).name = names{j};
            
            h3(1) = errorbar(xH, rH, rH_eb, 'marker', markers{1}, 'color', colorsHL{1}, 'linestyle', 'none');
            hold on
            h3(2) = errorbar(xL, rL, rL_eb, 'marker', markers{2}, 'color', colorsHL{2}, 'linestyle', 'none');
            
            
            %             if(nstim==2)
            %                 for i=1:length(rate)
            %                     n = nstim * (i-1);
            %                     for q=1:nstim
            %                         h3(q+n) = errorbar(xup{i}(:,q), rate{i}(:,q), reb{i}(:,q), 'marker', markers{i}, 'color', colors1{i}, 'linestyle', 'none');
            %                         hold on
            %                     end
            %                 end
            %             end
            title(names{j});
            ax = gca;
            ax.Box = 'off';
            ah(kk).handles = h3;
            
            if(nstim==2)
                ah(kk+1).pos = pos{kk+1}; %#ok<*AGROW>
                ah(kk+1).ax = axes('Position', pos{kk+1});
                ah(kk+1).name = names{j};
                h4(1) = errorbar(xHo, rHo, rHo_eb, 'marker', markers{1}, 'color', colorsBR{1}, 'linestyle', 'none');
                hold on
                h4(2) = errorbar(xLo, rLo, rLo_eb, 'marker', markers{1}, 'color', colorsBR{2}, 'linestyle', 'none');
                h4(3) = errorbar(xHl, rHl, rHl_eb, 'marker', markers{2}, 'color', colorsBR{3}, 'linestyle', 'none');
                h4(4) = errorbar(xLl, rLl, rLl_eb, 'marker', markers{1}, 'color', colorsBR{4}, 'linestyle', 'none');
                title(names{j});
                ax = gca;
                ax.Box = 'off';
                ah(kk+1).handles = h4;
            end
            
            %     set(ah(j).ax, 'Children', h);
            
            
        end
        
        
        
    case 'CO2'
        
        SF = 60;
        
        k = 0;
        for j = 1:2:length(btdstruct)
            
            k = k+1;
            ah(k).pos = pos{k}; %#ok<*AGROW>
            ah(k).ax = axes('Position', pos{k});
            ah(k).name = names{k};
            
            
            low = btdstruct(j).Rates.rf;
            high = btdstruct(j+1).Rates.rf;
            lx_low = low.lx;
            lx_high = high.lx;
            rate_low = SF*low.rate;
            rate_eb_low = SF*low.rate_eb;
            rate_high = SF*high.rate;
            rate_eb_high = SF*high.rate_eb;
            
            rateROG_low = SF*exp(polyval(low.rateFitQuad, lx_low));
            rateROG_high = SF*exp(polyval(high.rateFitQuad, lx_low));
            
            
            h(1) = errorbar (lx_low, rate_low, rate_eb_low, 'color', colors{k}{1}, 'marker', 'o', 'linestyle', 'none');
            hold all;
            h(2) = plot (lx_low, rateROG_low, 'color', 'k', 'linestyle', ':', 'linewidth', 1.5);
            
            h(3) = errorbar (lx_high, rate_high, rate_eb_high, 'color', colors{k}{2}, 'marker', 's', 'linestyle', 'none');
            h(4) = plot (lx_high, rateROG_high, 'color', 'k', 'linestyle', '--', 'linewidth', 1.5);
            
            title(names{k});
            
            ax = gca;
            ax.Box = 'off';
            ax.YLabel.String = 'turn rate (min^{-1})';
            ax.YLabel.FontSize = 10;
            ax.XLabel.String = 'Filtered Stimulus Value';
            ax.XLabel.FontSize = 10;
            
            %         set(ax, 'YScale', 'log');
            %         ax.YLim = [0 10];
            %         ax.YTick = [0 1 2 5 10];
            %         ax.YTickLabel = {'0', '1', '2', '5', '10'};
            
            ah(k).handles = h;
            
            
            %     set(ah(j).ax, 'Children', h);
            
            
        end
        
        
    case 'TRate'
        
        
        for j=1:length(btdstruct)
            ah(j).pos = pos{j}; %#ok<*AGROW>
            ah(j).ax = axes('Position', pos{j});
            ah(j).name = names{j};
            
            btd = btdstruct(j);
            tx = btd.histt(1:end-1);
            trate_on = btd.trate_on;
            trate_on_eb = btd.trate_on_eb;
            trate_off = btd.trate_off;
            trate_off_eb = btd.trate_off_eb;
            
            tr = trate_off;
            tr_eb = trate_off_eb;
            
            h = shadedErrorPlot (tx, tr, tr_eb, tr_eb, colors{j});
            
            hold (ah(j).ax, 'on');
            tm1 = max(tr(tx < max(tx) / 2));
            tm2 = max(tr(tx > max(tx) / 2));
            x1 = find(tr == tm1);
            x2 = find(tr == tm2);
            
            
            plot([tx(x1) tx(x2)], [tm1 tm1], 'k:', [tx(x1) tx(x2)], [tm2 tm2], 'k:');
            
            hold (ah(j).ax, 'off');
            
            ah(j).handles = h;
            title(names{j});
            
            ax = gca;
            ax.Box = 'off';
            
        end
        
        
    case 'blah'
        
        n = 1;
        for j=1:length(btdstruct)/2
            ah(j).pos = pos{j}; %#ok<*AGROW>
            ah(j).ax = axes('Position', pos{j});
            ah(j).name = names{j};
            
            btd_l = btdstruct(n);
            btd_h = btdstruct(n+1);
            lx_low = btd_l.lx{1};
            rate_low = btd_l.rate{1};
            rate_eb_low = btd_l.rate_eb{1};
            rateROG_low = btd_l.rateExp2{1};
            rateExp_low = btd_l.rateExp1{1};
            lx_high = btd_h.lx{1};
            rate_high = btd_h.rate{1};
            rate_eb_high = btd_h.rate_eb{1};
            rateROG_high = btd_h.rateExp2{1};
            rateExp_high = btd_h.rateExp1{1};
            h(1) = errorbar (lx_low, rate_low, rate_eb_low, 'color', colors{j}{1}, 'marker', 'o', 'linestyle', 'none');
            hold all;
            h(2) = plot (lx_low, rateROG_low, 'color', colors{j}{1}, 'linestyle', '--');
            
            h(3) = errorbar (lx_high, rate_high, rate_eb_high, 'color', colors{j}{2}, 'marker', 's', 'linestyle', 'none');
            h(4) = plot (lx_high, rateROG_high, 'color', colors{j}{2}, 'linestyle', '--');
            
            title(names{j});
            
            ax = gca;
            ax.Box = 'off';
            
            n = n+2;
            
            %         set(ax, 'YScale', 'log');
        end
        
end