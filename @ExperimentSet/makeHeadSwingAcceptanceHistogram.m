function [h,eb] = makeHeadSwingAcceptanceHistogram(eset, fieldname, fieldaxis, varargin)
% make a histogram of acceptance ratio vs. fieldname
% h = makeHeadSwingAcceptanceHistogram(eset, fieldname, fieldaxis, varargin)
%
% make a histogram of acceptance ratio vs. fieldname
% each headswing maps one value of fieldname (mean(fieldname(hs.inds))) to 
% a 0 (rejected) or 1 (accepted)
% we then form a graph of acceptance probability vs. fieldname within bins
% specified by fieldaxis
%
% to use value at start/end of headswing rather than mean value, pass
% 'atstart',true or 'atend', true as optional arguments;
%
% outputs:
%   H: (optional), the histogram;  if no output arguments, plots histogram
% inputs:
%   ESET: a member of the ExperimentSet class
%   FIELDNAME: the name of the field over which to make the acceptance
%   histogram
%   FIELDAXIS: determines the binning for the histogram
%   VARARGIN: optional parameter/value pairs
%       to use value at start/end of headswing rather than mean value, pass
%       'atstart',true or 'atend', true as optional arguments;

atstart = false;
atend = false;
varargin = assignApplicable(varargin);


t = [eset.expt.track];
hs = [t.headSwing];
acc = [hs.accepted];
%size(acc)
if (atstart)
    mfv = eset.gatherField(fieldname, 'hs', 'start');
else
    if (atend)
        mfv = eset.gatherField(fieldname, 'hs', 'end');
    else
        mfv = eset.gatherField(fieldname, 'hs','mean');
    end
end
%size(mfv)
n = hist (mfv, fieldaxis);
h1 = hist(mfv(logical(acc)), fieldaxis) ./ n;
eb = h1.*(1-h1)./sqrt(n);

if (nargout > 0)
    h = h1;
else
    errorbar (fieldaxis, h1, eb, varargin{:}); xlabel(fieldname); ylabel ('head swing acceptance');
    if atstart
        t = 'headswing acceptance vs. value at start of headswing'; %#ok<UNRCH>
    else if atend
            t = 'headwing acceptance vs. value at end of headswing'; %#ok<UNRCH>
        else
            t = 'headswing acceptance vs. mean value during headswing';
        end
    end
    title ([eset.defaultTitle ': ' t]);
end

