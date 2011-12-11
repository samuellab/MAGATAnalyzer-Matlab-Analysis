function calculateDerivedQuantity(expt, quantityNames, recalculate)
%calculates the derived quantity for each track using the same derivation rules
%function calculateDerivedQuantity(expt, quantityNames, recalculate)
%expt.calculateDerivedQuantity(quantityNames, recalculate)
%
%Assigns expt.dr to each track (overwriting any existing different
%derivation rules)
%then calls calculateDerivedQuantity on each track
%
%ouputs: none
%inputs:
%EXPT: a member of the experiment class
%QUANTITYNAMES: the names of quantities to calculate
%RECALCULATE: whether to recalculate the quantities if they are already
%calculated (optional, defaults to false)

existsAndDefault('recalculate', false);
[expt.track.dr] = deal(expt.dr);
expt.executeTrackFunction('calculateDerivedQuantity', quantityNames, recalculate);