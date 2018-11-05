
% This script can be used to plot different things from our Variance Adaptation paper

% Load the appropriate matfile and run the section that you need
% Each section calls VarAdaptPlots.m with different parameters and sets of experiments


%%
% run this first to set all the colors and figure positions

RedsHL = {1/255*[243 127 129], 1/255*[127 20 22]};
BluesHL = {1/255*[74 137 165], 1/255*[57 70 156]};
PurpleHL = {[.8 0.2 0.6], [.4 0 0.6]};
CyanHL = {[0 .2 .4], [0 .8 .8]};
Red = [1 0 0];
Blue = [0 0 1];
Cyan = [1 0 0.6];

RedAll = {[1 0 0], [1 .8 .8]};
BlueAll = {[0 0 1], [.8 .8 1]};
PurpleAll = {[.6 0 0.4], [.8 0.6 0.8]};
CyanAll = {[0 .4 .6], [.6 -8.2 .8]};

n = get(gcf);
if(n.Number == 1)
    fignum = 1;
else
    fignum = n.Number + -1;
end
[adim,po] = blank8x10Figure(fignum);
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


%%
[adim,po] = blank8x10Figure(fignum);

btds_40 = [berlin_40, or42a_40];
btds_120 = [berlin_120, or42a_120];

SF_120 = VarAdaptPlots(btds_120, {'', ''}, {{BlueAll}, {RedAll}}, {pos{3}, pos3{3}}, 'ScaleFactorVTime', 1, 'Bayes', '');
SF_40 = VarAdaptPlots(btds_40, {'', ''}, {{BlueAll}, {RedAll}}, {pos{2}, pos3{2}}, 'ScaleFactorVTime', 1, 'Bayes', '');

SF_120(1).ax.XLim = [-3 30];
SF_120(2).ax.XLim = [-3 30];
SF_40(1).ax.XLim = [-3 30];
SF_40(2).ax.XLim = [-3 30];

%%
% Unisensory, T=120s (Figure 1)

[adim,po] = blank8x10Figure(fignum);

log = 1;
fitdegree = 2; % 0 means don't plot the rate fits; 1,2 mean plot exponential fit with 1st or 2nd deg. polynomial

load('Z:\Var.Adapt Ruben\Matfiles\Uni-Sensory\T120s\berlin.mat'); % this is berlin_120
load('Z:\Var.Adapt Ruben\Matfiles\Uni-Sensory\T120s\or42a.mat'); % this is or42a_120


btds = [berlin_120, or42a_120];
lines_120 = {'berlin', 'or42a'};

tta = VarAdaptPlots(btds, lines_120, {BluesHL, RedsHL}, {pos{1}, pos3{1}}, 'TTA');
RF = VarAdaptPlots(btds, lines_120, {BluesHL, RedsHL}, {pos{2}, pos3{2}}, '1stimRF', 1, [], [], [], [], [], 1, fitdegree);
SF = VarAdaptPlots(btds, {'', ''}, {{BlueAll}, {RedAll}}, {pos{3}, pos3{3}}, 'ScaleFactorVTime', 1, 'Bayes', '');

set([RF(1).ax], 'YLim', [0 18]);
set([RF(2).ax], 'YLim', [0 12]);
set([RF(2).ax], 'YTick', [0 5 10]);

if(log)
    set([RF.ax], 'YScale', 'log');
    set([RF(1).ax], 'YLim', [0.5 18]);
    set([RF(2).ax], 'YLim', [0.5 12]);

    set([RF.ax], 'YTick', [0 1 2 5 10 15 20]);
    set([RF.ax], 'YTickLabel', {'0', '1', '2', '5', '10', '15', '20'});

end

%%
% UniSensory, T=40s (Figure 3)

n = get(gcf);
if(n.Number == 1)
    fignum = 1;
else
    fignum = n.Number + 1;
end
[adim,po] = blank8x10Figure(fignum);

log = 1;

% load matfiles from:
% Z:\Var.Adapt Ruben\Matfiles\Uni-Sensory\T40s\dt=0.1\eti\Quad Rate\

pdegree = 2;
btds = berlin_40;
colHL = RedsHL;
colAll = RedAll;

tta = VarAdaptPlots(btds, {''}, {colHL}, {pos2{1}}, 'TTA');
RF = VarAdaptPlots(btds, {''}, {colHL}, {pos2{2}}, '1stimRF', 0, [], [], [], [], [], 1, pdegree);
SF = VarAdaptPlots(btds, {''}, {{colAll}}, {pos3{1}}, 'ScaleFactorVTime', 0.1, 'Bayes', 'LowHigh');

if(log)
    set([RF.ax], 'YScale', 'log');
    set([RF.ax], 'YLim', [0.5 12]);
    set([RF.ax], 'YTick', [0 1 2 5 10 15 20]);
    set([RF.ax], 'YTickLabel', {'0', '1', '2', '5', '10', '15', '20'});

end




%% CO2 (Figure 7)

% load matfiles from:
% Z:\Var.Adapt Ruben\Matfiles\Uni-Sensory\Co2+Noise\Gr21a and
% Z:\Var.Adapt Ruben\Matfiles\Uni-Sensory\Co2+Noise\Or42a

log = 1;

co2_21a = VarAdaptPlots([gr21a_R2, gr21a_R8, gr21a_B2, gr21a_B8], {'', ''}, {RedsHL, BluesHL}, {pos{1}, pos{3}}, 'CO2');
co2_42a = VarAdaptPlots([or42a_R2, or42a_R8, or42a_B2, or42a_B8], {'', ''}, {RedsHL, BluesHL}, {pos2{1}, pos2{3}}, 'CO2');
% co2_42a(1).ax.XLabel.String = 'Filtered Stim. Value';
% co2_42a(2).ax.XLabel.String = 'Filtered Stim. Value';
co2_42a(1).ax.XAxis.Color = 'r';
co2_21a(1).ax.XAxis.Color = 'r';
co2_42a(2).ax.XAxis.Color = 'b';
co2_21a(2).ax.XAxis.Color = 'b';

% co2_21a(2).ax.YLim(2) = 25;

if(log)
    set([co2_21a.ax], 'YScale', 'log');
    set([co2_42a.ax], 'YScale', 'log');
    
    set(co2_21a(1).ax, 'YLim', [0.5 12]);
    set(co2_42a(1).ax, 'YLim', [0.5 12]);
    set(co2_21a(2).ax, 'YLim', [0.5 18]);
    set(co2_42a(2).ax, 'YLim', [0.5 18]);
    
    set([co2_21a.ax], 'YTick', [0 1 2 5 10 15 20]);
    set([co2_21a.ax], 'YTickLabel', {'0', '1', '2', '5', '10', '15', '20'});
    set([co2_42a.ax], 'YTick', [0 1 2 5 10 15 20]);
    set([co2_42a.ax], 'YTickLabel', {'0', '1', '2', '5', '10', '15', '20'});
end



%%
%Var Step + Sig. Step (not used in paper)
lines1 = {'berlin', 'berlin','or42a', 'or35a'};
lines2 = {'RN+Bstep64', 'BN+Rstep64', 'BN+Rstep96'};

ah_step1 = VarAdaptPlots([berlin_step, berlin_step2, or42a_step, or35a_step], lines1, {Blue, Blue, Red, Red}, {pos{1}, pos{2}, pos2{1}, pos2{2}}, 'TRate');

% ah_step2 = VarStepRatePlots([RN_Bstep64, BN_Rstep64, BN_Rstep96], lines2, {Blue, Red, Red}, {pos3{1}, pos3{2}, pos3{3}}, 'TRate');


%%
% Var Ramps (Figure 4)

[adim,po] = blank8x10Figure();
log = 1;
deltaT = 1;
pdegree = 2;

% load matfiles from:
% Z:\Var.Adapt Ruben\Matfiles\Uni-Sensory\Ramps\dt=1\Quad Rate

lines = {''};
btds = [or42a_ramp];
col = 'Red';
colHL = RedsHL;

SF = VarAdaptPlots(btds, lines, {col}, {pos{1}, pos{2}}, 'VarRamp', deltaT, 'Bayes');

% RF = VarAdaptPlots(btds, lines, {colHL}, {pos{3}, pos3{3}}, '1stimScaledRF', deltaT, [], [], [], [], [], [], pdegree);

ah_ramp(1).ax.XLabel.String = 'time in cycle (s)';

if(log)
    set([RF.ax], 'YScale', 'log');
    set([RF.ax], 'YLim', [0.5 12]);

    set([RF.ax], 'YTick', [0 1 2 5 10 15 20]);
    set([RF.ax], 'YTickLabel', {'0', '1', '2', '5', '10', '15', '20'});

end

%%
% Corr Ramp (not used in paper)
% figure()

[adim,po] = blank8x10Figure(fignum);

log = 1;
deltaT = 1;
lines = {''};
btds = rbcorr_ramp;

col = PurpleAll;

SF = VarAdaptPlots(btds, lines, {col}, {pos{1}, pos{2}}, 'VarRamp', deltaT, 'Bayes');
RF = VarAdaptPlots(btds, {''}, {PurpleHL}, {pos3{1}}, '1stimRF', [], '', [], [], [], [], 1);
RF2 = VarAdaptPlots(btds, {''}, {PurpleHL}, {pos3{2}}, '1stimRF', [], '', [], [], [], [], 2);
RF.ax.XAxis.Color = PurpleAll{1};

if(log)
    set([RF.ax], 'YScale', 'log'); 
    set([RF.ax], 'YLim', [0.5 12]);
    set([RF2.ax], 'YScale', 'log');
    set([RF2.ax], 'YLim', [0.5 12]);
    
    set([RF.ax], 'YTick', [0 1 2 5 10 15 20]);
    set([RF.ax], 'YTickLabel', {'0', '1', '2', '5', '10', '15', '20'});
    set([RF2.ax], 'YTick', [0 1 2 5 10 15 20]);
    set([RF2.ax], 'YTickLabel', {'0', '1', '2', '5', '10', '15', '20'});

end

%%
% Corr Step (Figure 6, bottom)

[adim,po] = blank8x10Figure(fignum);
log = 1;

PurpleAll = {[.6 0 0.4], [.8 0.6 0.8]};
colHL = PurpleHL;

% load matfile from:
% load('Z:\Var.Adapt Ruben\Matfiles\Multi-Sensory\Correlated\rbcorr_Quad.mat');

deltaT = 0.1;
btd = rbcorr;

SF = VarAdaptPlots(btd, {''}, {{PurpleAll PurpleAll}}, {pos{1}}, 'ScaleFactorVTime', deltaT, 'Bayes', [], [], [], [] , 1);
RF = VarAdaptPlots(btd, {''}, {PurpleHL}, {pos{2}}, '1stimRF', [], '', [], [], [], [], 1, 1);

RF.ax.XAxis.Color = PurpleAll{1};
RF2.ax.XAxis.Color = PurpleAll{1};

if(log)
    set([RF.ax], 'YScale', 'log');
    set([RF.ax], 'YLim', [0.5 12]);
        
    set([RF2.ax], 'YScale', 'log');
%     set([RF2.ax], 'YScale', 'log');
    
    set([RF.ax], 'YTick', [0 1 2 5 10 15 20]);
    set([RF.ax], 'YTickLabel', {'0', '1', '2', '5', '10', '15', '20'});
    set([RF2.ax], 'YTick', [0 1 2 5 10 20]);
    set([RF2.ax], 'YTickLabel', {'0', '1', '2', '5', '10', '20'});

end

%%
% Multi-Sensory Variance Switching, T=40s (Figure 6, top 2 rows)

[adim,po] = blank8x10Figure();
log = 1;


% matfiles in: Z:\Var.Adapt Ruben\Matfiles\Multi-Sensory\VarStep_T40\dt=0.1\Quad Rate

rcbsw = {'', '', ''};
% bcrsw = {'Blue s=1', 'Blue s=2', 'Blue s=3'};
% btds = [r1b13 r2b13 r3b13];
btds = [b2r13 r2b13];
pdegree = 1;

SF = VarAdaptPlots(btds, rcbsw, {RedAll{1} BlueAll{1}}, {pos{1}, pos2{1}, pos3{1}},  '2stimPPF', 0.1, 'Bayes', [],[],[],[],[],[],'');
RF = VarAdaptPlots(btds, rcbsw, {RedsHL, BluesHL}, {{pos{3}, pos{2}}, {pos2{3} pos2{2}}, {pos3{3} pos3{2}}}, '2stimRF', [], [], [], [], [], [], [], pdegree);


RF{1}(1).ax.YLim = [0 12];
RF{2}(1).ax.YLim = [0 18];
RF{1}(2).ax.YLim = [0 12];
RF{2}(2).ax.YLim = [0 18];
% RF{1}(3).ax.YLim = [0 13];
% RF{2}(3).ax.YLim = [0 13];

RF{2}(1).ax.XLim = RF{2}(2).ax.XLim;
RF{1}(1).ax.XLim = RF{1}(2).ax.XLim;

if(log)

    set([RF{1}.ax], 'YScale', 'log');
    set([RF{2}.ax], 'YScale', 'log');

    set([RF{1}.ax], 'YTick', [0 1 2 5 10 15 20]);
    set([RF{1}.ax], 'YTickLabel', {'0', '1', '2', '5', '10', '15', '20'});
    set([RF{2}.ax], 'YTick', [0 1 2 5 10 15 20]);
    set([RF{2}.ax], 'YTickLabel', {'0', '1', '2', '5', '10', '15', '20'});

end

%%
% Bayes Optimal Estimates (Figure 5, b)

% this plots optimal estimates of the variance of the stimulus in or42a_40 or berlin_40, converted to alpha(t) vs. cycle time

% conversion from sigma to alpha: alpha(sigma) = 1/sqrt(sigma^2+sigma_0^2);
% sigma_0 = params as set below (found separately for or42a_40 and berlin_40)

% you need to add the btds before the Bayesian estimator to find the stimulus variance
% for example:
% load('Z:\Var.Adapt Ruben\Matfiles\Uni-Sensory\T40s\dt=0.1\eti\QuadRate\or42a.mat');
% or42a_40_data = BehaviorTriggeredData.loadBTDDirectory('Z:\Var.Adapt Ruben\BTDs\By Expt Type\Var. Step T=40\Or42a');
% or42a_40.btd = or42a_40_data.btd;

% now you can run this section

[adim,po] = blank8x10Figure();

deltaT = 0.1;


dts1 = [.2 .4 .7 1.5];
dts2 = [.2 .4 .85 1.5];

taus = [1 1 1 1 1];
taus1 = 12*taus;
taus2 = 6*taus;

params = {1.9};
params2 = {1.3};
col = RedAll;
col2 = BlueAll;
 
ah = VarAdaptPlots([or42a_40], {''}, {{col}}, {pos{2}}, 'ScaleFactorVTime', deltaT, 'Bayes', 'Optimal', dts1, taus1, params, []);
ah2 = VarAdaptPlots([berlin_40], {''}, {{col2}}, {pos{1}}, 'ScaleFactorVTime', deltaT, 'Bayes', 'Optimal', dts2, taus2, params2, []);







