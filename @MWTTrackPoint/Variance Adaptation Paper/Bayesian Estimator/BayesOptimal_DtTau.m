function SF = BayesOptimal_DtTau(btdstruct, dts, taus, isconvolved)
%function SF = BayesOptimal_DtTau(btdstruct, dts, taus, isconvolved)

% btdstruct needs to contain the btd field which contains the stimulus
% isconvolved = 0 or 1. whether the filtered or raw stimulus is used to estimate the variance.
% dts,taus: arrays of desired values


SF.dts = dts;
SF.taus = taus;


VarEst = @(dt, tau) AlphaEstimate(btdstruct, btdstruct.varops, 'eti', isconvolved, dt/tau, dt, 60:dt:1200, btdstruct.var.period, btdstruct.var.tshift);

for i=1:length(dts)
    for k=1:length(taus)
        tic
        SF.sim{i}{k} = VarEst(dts(i), taus(k));
        toc
    end
end