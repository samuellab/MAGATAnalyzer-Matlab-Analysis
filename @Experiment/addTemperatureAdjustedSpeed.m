function addTemperatureAdjustedSpeed (expt, varargin)
% function addTemperatureAdjustedSpeed (expt, varargin)
%
% calculates mean speed vs. temp, then creates adjusted speed to take out linear contribution
% see addAdjustedField for implementation details
% EXPT < expt

fn = {expt.globalQuantity.fieldname};
ind = strcmpi('temperature', fn);
if(any(ind))
    fn = fn{find(ind, 1, 'first')};
    expt.addAdjustedField(fn, 'speed');
    return;
end
ind = strcmpi('temp', fn);
if(any(ind))
    fn = fn{find(ind, 1, 'first')};
    expt.addAdjustedField(fn, 'speed');
    return;
end    