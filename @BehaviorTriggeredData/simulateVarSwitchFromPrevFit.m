function problemDescription = simulateVarSwitchFromPrevFit (btdstruct, osfit, oldpd, alphaLow, alphaHigh, nlarvae)
% problemDescription = simulateVarSwitchFromPrevFit (btdstruct, osfit, oldpd, alphaLow, alphaHigh, nlarvae)
%

btd = btdstruct.btd(min(oldpd.exprange):max(oldpd.exprange));
opstruct = btdstruct.varops;
opstruct.taxis = min(oldpd.runEti):oldpd.deltaT:max(oldpd.runEti);
opstruct.alphaLow = alphaLow;
opstruct.alphaHigh = alphaHigh;

problemDescription =  simulateVarSwitch (btd, opstruct, oldpd, oldpd.timeField, oldpd.timeRange, osfit.staticParams, nlarvae);

problemDescription.params_0 = osfit.staticParams;
problemDescription.alpha_0 = osfit.alpha;
problemDescription.tx = osfit.tx;
% problemDescription.tparams_0 = osfit.temporalParams;
fn = {'pad', 'Q_alpha', 'v0', 'w_0', 'maxreps', 'separateExperiments', 'period', 'tshift'};
for j = 1:length(fn)
    problemDescription.(fn{j}) = oldpd.(fn{j});
end