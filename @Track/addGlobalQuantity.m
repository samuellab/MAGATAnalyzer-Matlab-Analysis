function addGlobalQuantity (track, quantityName, time, value) 
%function addGlobalQuantity (track, quantityName, time, value) 
%
%this function will eventually migrate to
%addGlobalQuantity(track,globalQuantity) where globalQuantity is a member
%of the GlobalQuantity class
%
%adds a field called quantityName to track.dq
%value of this field is interp1(time, value, track.dq.eti, 'linear', NaN)
disp('Track.addGlobalQuantity is deprecated and will be removed;  use addQuantityToTrack (gq, track) in GlobalQuantity class instead');

track.dq.(quantityName) = interp1(time, value, track.getDerivedQuantity('eti'), 'linear', NaN);

