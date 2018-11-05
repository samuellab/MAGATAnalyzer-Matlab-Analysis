function [ RF ] = plotMultiSensoryFits( fitstruct, varstruct, opstruct, color)


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

positions = {pos{1}, pos2{1}, pos2{2}, pos2{3}, pos3{1}, pos3{2}, pos3{3}};

lograte = 1;

for k=1:length(fitstruct.ratefuns)
    
    for i=1:length(fitstruct.nLL_meanR)
        if(strncmpi(fitstruct.models{k}, 'ind', 3))
            ao = fitstruct.rP{k}{i}(1);
            bo = fitstruct.rP{k}{i}(2);
            co = fitstruct.rP{k}{i}(3);
            al = fitstruct.rP{k}{i}(4);
            bl = fitstruct.rP{k}{i}(5);
            cl = fitstruct.rP{k}{i}(6);
        elseif(strncmpi(fitstruct.models{k}, 'mult', 3))
            a(i) = fitstruct.rP{k}{i}(1);
            b(i) = fitstruct.rP{k}{i}(2);
            c(i) = fitstruct.rP{k}{i}(3);
            d(i) = fitstruct.rP{k}{i}(4);
            e(i) = fitstruct.rP{k}{i}(5);
        elseif(strncmpi(fitstruct.models{k}, 'add', 3))
            theta(i) = fitstruct.rP{k}{i}(1);
            add_quad(i) = fitstruct.rP{k}{i}(2);
            add_lin(i) = fitstruct.rP{k}{i}(3);
            r0(i) = fitstruct.rP{k}{i}(4);
        end
    end
    
end

model_params = {};
if(any(strncmpi(fitstruct.models, 'ind', 3)))
    params_ind = [mean(ao) mean(bo) mean(co) mean(al) mean(bl) mean(cl)];
    model_params{end+1} = params_ind;
end
if(any(strncmpi(fitstruct.models, 'add', 3)))
    params_add = [mean(theta) mean(add_quad) mean(add_lin) mean(r0)];
    model_params{end+1} = params_add;
end
if(any(strncmpi(fitstruct.models, 'mult', 3)))
    params_mult = [mean(a) mean(b) mean(c) mean(d) mean(e)];
    model_params{end+1} = params_mult;
end


nBins = opstruct.numLxBins;
deltaT = median(diff(varstruct(1).fullensemble.eti));

for i=1:length(varstruct)
    td(:,i) = varstruct(i).turn.x_conv;
    rd(:,i) = varstruct(i).noturn.x_conv;
end

Rates2D = turnRates2D(td, rd, deltaT, nBins, fitstruct.ratefuns, model_params);

colormap jet
paxis = [0 .2];
pticks = [0 .1 .2];

RF{1}.pos = positions{1}; %#ok<*AGROW>
RF{1}.ax = axes('Position', positions{1});
pcolor(Rates2D.lx(:,1),Rates2D.lx(:,2), Rates2D.thist./(.2*.2*sum(sum(Rates2D.thist)))); axis square; shading interp; colorbar
RF{1}.ax.YAxis.TickValues = [-2 0 2];
clear caxis
caxis(gca, paxis); a = colorbar; a.Ticks = pticks;

nm = length(model_params);
for k=1:nm
    RF{k+1}.pos = positions{k+1};
    RF{k+1}.ax = axes('Position', positions{k+1});
    pcolor(Rates2D.lx(:,1),Rates2D.lx(:,2), Rates2D.pred_thist_ROG{k}); axis square; shading interp;
    title(fitstruct.models{k});
    clear caxis
    caxis(gca, paxis); a = colorbar; a.Ticks = pticks;
end
set(gcf, 'renderer', 'painters');
return

%     RF{k+nm+1}.ax.XAxis.TickValues = [-2 0 2];
%     RF{k+nm+1}.ax.YAxis.TickValues = [-2 0 2];
%     RF{k+nm+1}.pos = positions{k+nm+1};
%     RF{k+nm+1}.ax = axes('Position', positions{k+nm+1});
%     pcolor(Rates2D.lx(:,1),Rates2D.lx(:,2), Rates2D.pred_thist{k}); axis square; shading interp;
%     clear caxis
%     caxis(gca, paxis); a = colorbar; a.Ticks = pticks;
%     RF{k+nm+1}.ax.XAxis.TickValues = [-2 0 2];
%     RF{k+nm+1}.ax.YAxis.TickValues = [-2 0 2];
% end

% set(gcf, 'renderer', 'opengl');

% return

end

function rf = turnRates2D(td, rd, deltaT, nBins, ratefuns, params)

lx(:,1) = linspace(-3*std(rd(:,1)), 3*std(rd(:,1)), nBins);
lx(:,2) = linspace(-3*std(rd(:,2)), 3*std(rd(:,2)), nBins);
% lx(:,1) = linspace(-2.5, 2.5, nBins);
% lx(:,2) = linspace(-2.5, 2.5, nBins);

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

rf.muO_turn = mean(td(:,1));
rf.muL_turn = mean(td(:,2));
rf.muOL_turn = mean(td(:,1).*td(:,2),1);

xyall_gaussian = reshape(mvnpdf(lxdata, mean(rd,1), [std(rd(:,1)) 0;0 std(rd(:,2))]), size(lxx));
% xyall_gaussian = reshape(mvnpdf(lxdata, [0 0], eye(2)), size(lxx));

if(~isempty(ratefuns))
    for k=1:length(ratefuns)
        
        pred_rate{k} = reshape(ratefuns{k}(params{k}, lxdata), size(lxx));
        

        pred_thist_ROG{k} = pred_rate{k}.*xyall_gaussian;
        pred_thist{k} = pred_rate{k}.*rhist;
        
        rf.pred_thist_ROG{k} = pred_thist_ROG{k}./(prod(median(diff(lx,1)))*sum(sum(pred_thist_ROG{k})));
        rf.pred_thist{k} = pred_thist{k}./(prod(median(diff(lx,1)))*sum(sum(pred_thist{k})));
        
        rf.muO_turn_pred(k) = (prod(median(diff(lx,1))))*sum(sum(repmat(lx(:,1)',[length(lx(:,2)),1]).*rf.pred_thist{k}));
        rf.muL_turn_pred(k) = (prod(median(diff(lx,1))))*sum(sum( repmat(lx(:,2),[1,length(lx(:,1))]).*rf.pred_thist{k}));
        rf.muOL_turn_pred(k) = (prod(median(diff(lx,1))))*sum(sum(repmat(lx(:,1)',[length(lx(:,2)),1]).*repmat(lx(:,2), [1, length(lx(:,1))]).*rf.pred_thist{k}));
        
        rf.muO_turn_pred_ROG(k) = (prod(median(diff(lx,1))))*sum(sum(repmat(lx(:,1)',[length(lx(:,2)),1]).*rf.pred_thist_ROG{k}));
        rf.muL_turn_pred_ROG(k) = (prod(median(diff(lx,1))))*sum(sum( repmat(lx(:,2),[1,length(lx(:,1))]).*rf.pred_thist_ROG{k}));
        rf.muOL_turn_pred_ROG(k) = (prod(median(diff(lx,1))))*sum(sum(repmat(lx(:,1)',[length(lx(:,2)),1]).*repmat(lx(:,2), [1, length(lx(:,1))]).*rf.pred_thist_ROG{k}));
    end
end


end
    
    
    
    

