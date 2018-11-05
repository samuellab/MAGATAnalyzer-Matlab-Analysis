function addTemperatureAdjustedSpeed (eset, varargin)
% function addTemperatureAdjustedSpeed (expt, varargin)
%
% calculates mean speed vs. temp, then creates adjusted speed to take out linear contribution
% see addAdjustedField for implementation details
% EXPT < expt

fn = {eset.expt(1).globalQuantity.fieldname};
ind = strcmpi('temperature', fn);
if(any(ind))
    fn = fn{find(ind, 1, 'first')};
    eset.addAdjustedField(fn, 'speed');
    return;
end
ind = strcmpi('temp', fn);
if(any(ind))
    fn = fn{find(ind, 1, 'first')};
    eset.addAdjustedField(fn, 'speed');
    return;
end    