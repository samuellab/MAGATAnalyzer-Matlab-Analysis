function addStandardizedField(expt, fieldname, varargin)
% creates a standardized field for a given quantity
% function addStandardizedField(expt, fieldname, varargin)
%
% if x(j) is the value of the field then 
% standardized_x(j) = (x(j) - <x>) / std(x), where the mean and standard
% deviation are taken across the whole experiment
if (length(expt) > 1)
    for j = 1:length(expt)
        expt(j).addStandardizedField(fieldname, varargin{:});
    end
    return;
end

if (iscell(fieldname))
    for j = 1:length(fieldname)
        expt.addStandardizedField(fieldname{j}, varargin{:});
    end
    return;
end

newfieldname = ['standardized_' fieldname];
qv = (expt.gatherField(fieldname));
meanf = mean(qv);
stdf = std(qv);

gq = GlobalQuantity();
gq.xField = fieldname;
gq.fieldname = newfieldname;
gq.xData = meanf;
gq.yData = stdf;
gq.derivationMethod = @(xin, xData, yData) (xin - xData)./yData;

expt.addGlobalQuantity(gq);

end

