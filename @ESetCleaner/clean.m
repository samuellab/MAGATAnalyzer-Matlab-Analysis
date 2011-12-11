function clean(ecl, eset)
% function clean(ecl, eset)
%
% removes tracks that likely contain invalid data, according to criteria 
% stored in ecl
%
% outputs: none
% inputs:
%   ecl < ESetCleaner
%   eset < ExperimentSet

if (ecl.askFirst)
    if (ecl.getReport(eset) == 0)
        %nothing to clean
        return;
    end
    
    key = input ('continue with cleaning? : y/[n]', 's');
    if (isempty(key) || lower(key(1)) ~= 'y')
        return;
    end
end

for j = 1:length(eset.expt)

    sp = eset.expt(j).gatherField('speed', 'mean');

    if (eset.expt(j).track(1).validDQName('ihtValid'))
        htv = eset.expt(j).gatherField('ihtValid', 'mean');
    else
        htv = ones(size(sp));
    end
    dst = eset.expt(j).evaluateTrackExpression('max(sqrt(sum(track(1).getDerivedQuantity(''displacement'').^2)))');
    npts = eset.expt(j).gatherField('npts');
    
    pa = eset.expt(j).evaluateTrackExpression('max(abs(unwrap(track.getDerivedQuantity(''theta''))))');
    et = eset.expt(j).evaluateTrackExpression('max(track.getDerivedQuantity(''eti'')) - min(track.getDerivedQuantity(''eti''))');
   % pa = pa / (2*pi);
    %et = et / 60;
    rpm = pa./et * 60 / (2*pi);
    
    valid = htv >= ecl.minHTValid & npts >= ecl.minPts & dst >= ecl.minDist & sp >= ecl.minSpeed & (rpm < ecl.rpmCut | pa < ecl.minRevCut * 2*pi);
    eset.expt(j).track = eset.expt(j).track(valid);
end