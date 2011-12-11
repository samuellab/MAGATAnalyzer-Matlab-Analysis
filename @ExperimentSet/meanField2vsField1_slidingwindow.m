function [x,meany,standarderror,standarddeviation] = meanField2vsField1_slidingwindow (eset, field1name, field2name, field1centers, field1binsize, windowType, varargin)
% the mean value of points in field2 corresponding to values of field1
% function [x,meany,standarderror,standarddeviation] = meanField2vsField1
% (eset, field1name, field2name, field1axis, varargin)
%
% gathers two fields (with varargin arguments for both) and calls meanyvsx on them
% if no return arguments are requested, plots the result
% note that field1axis is bin edges [edge1,edge2) and x is the mean value
% of field1 in that bin
%
% outputs:
% X: the x-axis, the mean value of field1 for all points in each bin
% MEANY: the mean value of field2 for all points in each bin
% STANDARDERROR: the standard deviation of field2 in each bin / 
%   sqrt of number of points in each bin
% STANDARDDEVIATION: the standard deviation of field2 in each bin
% inputs:
% ESET: a member of the ExperimentSet class
% FIELD1NAME: the name of the independent field (must be 1D)
% FIELD2NAME: the name of the dependent field (may be ND)
% FIELD1AXIS: bin edges [edge1,edge2); thus length(x) =
%       length(field1axis) - 1; and x(j) is the mean value of all points
%       with field1 >= field1axis(j) and field1 < field1axis(j+1)
% VARARGIN:
%       any option that can be passed to ExperimentSet/gatherField
%       'inds', inds will take value only from inds in this
%       fashion (same for ydata/field2name):
%       xdata = eset.gatherField(field1name, varargin{:});
%       xdata = xdata(inds)
%       'field1dim' - if field1 is ND instead of 1D, choose this dimension
%                   - e.g. if field1name is 'loc', pass 'field1dim',1 for x
%                     position
%       'polar', true/false -- adjusts xdata for polar histogram
%       'autocorr_tau' - autocorrelation time constant (in seconds)
%                      <q(t0),q(t+t0)>_t0 ~= exp(-t/tau)
%       'timerange' -- only consider within this time period
inds = [];
field1dim = [];
autocorr_tau = -1;
timerange = [];
polar = false;
varargin = assignApplicable(varargin);
if (autocorr_tau < 0)
    if (eset.autocorr_tau <= 0)
        disp ('calculating and setting autocorrelation time constant');
        eset.setAutocorrTau;
    end
    autocorr_tau = eset.autocorr_tau;
end
xdata = eset.gatherField(field1name, varargin{:});
if ~isempty(field1dim)
    xdata = xdata(field1dim,:);
end

ydata = eset.gatherField(field2name, varargin{:});
if (~isempty(timerange))
    eti = eset.gatherField('eti', varargin{:});
    timinds = eti >= min(timerange) & eti <= max(timerange);
    if (~isempty(inds))
        if (max(inds) > 1)
            inds2 = false(size(timinds));
            inds2(inds) = true;
            inds = inds2;
        end
        inds = inds & timinds;
    else
        inds = timinds;
    end
end
if (~isempty(inds))
    xdata = xdata(inds); %xdata must be 1-D
    ydata = ydata(:,inds);
end
[x,meany,standarderror,standarddeviation] = meanyvsx_slidingwindow(xdata, ydata, field1centers, field1binsize, windowType, polar);
dti = eset.expt(1).track(1).dr.interpTime;
k = sqrt((1 + exp(-dti/autocorr_tau))/(1 - exp(-dti/autocorr_tau)));
standarderror = k*standarderror;
if (nargout == 0)
    shadedErrorPlot (x,meany,standarderror); title (eset.defaultTitle);
    xlabel(field1name);
    ylabel(['$<$' field2name '$>$']);
    embiggen();
end