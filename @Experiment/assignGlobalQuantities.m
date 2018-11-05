function assignGlobalQuantities(expt,varargin)
%reassigns all tracked global quantities to all tracks
%function assignGlobalQuantities(expt,varargin)
%expt.assignGlobalQuantities
%
%reassigns all tracked global quantities (see expt.globalQuantity)
%to all tracks
if (length(expt) > 1)
    for k = 1:length(expt)
        expt(k).assignGlobalQuantities;
    end
    return;
end

for k = 1:length(expt.globalQuantity);
    gq = expt.globalQuantity(k);
    for j = 1:length(expt.track)
        gq.addQuantityToTrack(expt.track(j));
    end
end

%{
for k = 1:length(expt.globalLookupTable);
    gq = expt.globalLookupTable(k);
    for j = 1:length(expt.track)
        gq.addQuantityToTrack(expt.track(j));
    end
end
%}