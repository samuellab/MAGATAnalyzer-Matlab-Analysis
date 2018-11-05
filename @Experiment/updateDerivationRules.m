function updateDerivationRules(expt, varargin)
% changes all derivation rules be the same as experiment's
% does not recalculate derived quantities
% function updateDerivationRules(expt)
% function updateDerivationRules(expt, dr)
% EXPT < experiment
% DR < DerivationRules (optional) - sets expt.dr = dr

if (length(expt) > 1)
    for j = 1:length(expt)
        expt(j).updateDerivationRules(varargin{:});
    end
    return;
end
if (length(varargin) >= 1)
    dr = varargin{1};
    if (isa(dr, 'DerivationRules'))
        expt.dr = dr;
    end
end
[expt.track.dr] = deal(expt.dr);