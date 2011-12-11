function [h,eb] = makeHistogram(eset, fieldname, fieldaxis, varargin)
% generates a histogram of expt.track.getDerivedQuantity(fieldname,varargin{:})
% function [h,eb] = makeHistogram(eset, fieldname, fieldaxis, varargin)
%
% generates a histogram of all values in expt.track.getDerivedQuantity(fieldname,varargin{:})
% if no arguments are specified, generates a plot of that histogram
%
% 
% outputs:
%   H: (optional), the histogram;  if no output arguments, plots histogram
% inputs:
%   ESET: a member of the ExperimentSet class
%   FIELDNAME: the name of the field over which to make the histogram 
%   FIELDAXIS: bin centers for the histogram
%   VARARGIN: 
%      anything to ExperimentSet.gatherField
%      optional parameter/value pairs
%      'r2d',true - when plotting, plot fieldaxis in degrees instead of radians
%      'polar',true - adjust data so that the 0 degree bin is the same size as the
%           next bin;  the last bin may be screwed by this, so be careful
%           in other words, passing polar, fieldaxis = deg2rad(0:30:330) is good
%           passing polar,true deg2rad(0:30:360) will nuke your last bin
%           it's necessary to pass polar as an argument if you are using
%           time chunks to calculate proper error bars and the data is
%           polar;  otherwise the averaging due to time chunks can cause
%           problems if the chunk passes through the x-axis
%      'autocorr_tau' - autocorrelation time constant (in seconds)
%                      <q(t0),q(t+t0)>_t0 ~= exp(-t/tau)

r2d = false;
polar = false;
%timechunk = [];
autocorr_tau = -1;
varargin = assignApplicable(varargin);
data = eset.gatherField(fieldname, varargin{:});
if (autocorr_tau < 0)
    if (eset.autocorr_tau <= 0)
        disp ('calculating and setting autocorrelation time constant');
        eset.setAutocorrTau;
    end
    autocorr_tau = eset.autocorr_tau;
end

if (polar)
    data = unwrap(data);
end
%{
if (~isempty(timechunk))
    t = eset.gatherField('eti', varargin{:});
    tx = min(t):timechunk:(max(t)+timechunk); %#ok<BDSCI>
    binsize = timechunk/eset.expt(1).dr.interpTime;
    [~,my,~,~,sumy] = meanyvsx(t,data,tx);
    nelems = sumy./my;
    howmany = nelems./binsize;
    data = my(isfinite(howmany) & howmany > 0);
    howmany = howmany(isfinite(howmany) & howmany > 0);
else
    howmany = ones(size(data));
end
%}

if (polar)
    c1 = fieldaxis(1);
    dt = fieldaxis(2) - c1;
    edge = c1 - dt/2;
    data = mod(data-edge,2*pi)+edge;
end
h1 = hist(data(isfinite(data)), fieldaxis);
n = sum(h1);
dti = eset.expt(1).track(1).dr.interpTime;
k = sqrt((1 + exp(-dti/autocorr_tau))/(1 - exp(-dti/autocorr_tau)));
eb = k*sqrt((h1/n).*(1-h1/n))*sqrt(n); %sigma = sqrt(Npq)
%{
low = 2*fieldaxis(1)-fieldaxis(2);
high = 2*fieldaxis(end)-fieldaxis(end-1);
fieldedges = 1/2 * ([low fieldaxis] + [fieldaxis high]);
[~,bin] = histc(data,fieldedges);
h1 = zeros(size(fieldaxis));

for j = 1:length(fieldaxis)
    h1(j) = sum(howmany(bin == j));
end
%}
if (nargout == 0)
    if (r2d)
        bar (rad2deg(fieldaxis), h1); hold on; title (eset.defaultTitle);
        errorbar (rad2deg(fieldaxis),  h1, eb,'k.', 'LineWidth', 3); hold off
    else
        bar (fieldaxis, h1); hold on; title (eset.defaultTitle);  
        errorbar (fieldaxis, h1, eb,'k.', 'LineWidth', 3); hold off
    end
    
    xlabel(fieldname); ylabel('counts'); embiggen();
else
    h = h1;
    
end
