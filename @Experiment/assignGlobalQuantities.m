function assignGlobalQuantities(expt,varargin)
%reassigns all tracked global quantities to all tracks
%function assignGlobalQuantities(expt,varargin)
%expt.assignGlobalQuantities
%
%reassigns all tracked global quantities (see expt.globalQuantity)
%to all tracks

for k = 1:length(expt.globalQuantity);
    gq = expt.globalQuantity(k);
    for j = 1:length(expt.track)
        gq.addQuantityToTrack(expt.track(j));
    end
end
