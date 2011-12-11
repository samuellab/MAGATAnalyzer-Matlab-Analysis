function setAutocorrTau (eset, varargin)
%function setAutocorrTau (eset, varargin)
%
%sets the autocorrelation time constant
%uses the autocorrelation of the field (default vnorm)
%and sets tau to be the first time s.t. ac(tau) <= 1/e
%
%optional arguments:
%'fieldname', fieldname (default 'theta') -- but note default call with no
%     additional arguments is different
%calls to ExperimentSet.autocorrelate, including
%'row', row number(s)
%if row = m, just finds the correlation of track.dq.(fieldnames)(m,:)\
%'inRuns', true/false
%if inRuns is true, we take the autocorrelation over the whole track, but
%first interpolate the fields over only the run indices
%this is useful if the fields are ill-defined between runs (e.g. velocity
%direction)
%'withinRuns', true/false
%if withinRuns is true, we find the correlation only within each run
%'isangle', true/false
%   if isangle is true, then we compute the correlation as cos (theta1 -
%   theta2) instead of as theta1*theta2
%
%default autocorrelation is 'fieldname', theta, 'isangle', true, 'inRuns',
%true
%and is performed if no optional args are passed in

eset.autocorr_tau = eset.getAutocorrTau(varargin{:});

    
