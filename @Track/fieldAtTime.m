function qvec = fieldAtTime(track, quantityName, elapsedTime) 
% gets the value(s) of track.dq.quantityName interpolated at elapsedTime(s)
% function qvec = fieldAtTime(track, quantityName, elapsedTime) 
%
% gets the value(s) of track.dq.quantityName interpolated at elapsedTime(s)
% NaN is returned for any times before or after the end of the track
%
% outputs: 
% QVEC: a kxN array of values
% inputs:
% TRACK: a member of the track class
% QUANTITYNAME: the name of the quantity (see track.getDerivedQuantity)
% ELAPSEDTIME: the times at which to get the values

val = transpose(track.getDerivedQuantity(quantityName));
x = transpose(track.getDerivedQuantity('eti'));

qvec = transpose(interp1(x, val, elapsedTime, 'linear', NaN));
