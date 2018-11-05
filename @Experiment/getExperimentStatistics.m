function es = getExperimentStatistics(expt, varargin)
if (length(expt) > 1)
    for j = 1:length(eset)
        es(j) = getExperimentStatistics(expt(j), varargin{:}); %#ok<AGROW>
    end
    return;
end


eti = expt.gatherField('eti');
timerange = [min(eti)-1 max(eti)+1];
es.timerange = timerange;

it = expt.gatherSubField('dr','interpTime');
dt = median(it);
if (any(it ~= dt))
    disp (['warning:  eset does not have homogenous interpolation times, instead range from ' num2str(min(it)) ' to ' num2str(max(it))]);
end


eet = expt.gatherField('eti');

es.numAnimals = max(histc(eet, min(eet):(100*dt):max(eet)))/100;
tx = min(eet):30:max(eet);
h = hist(eet, tx);
[~,I] = max(h);


sl  = expt.gatherField('spineLength', 'mean');
st = expt.evaluateTrackExpression('min(track.getDerivedQuantity(''eti''))');
et = expt.evaluateTrackExpression('max(track.getDerivedQuantity(''eti''))');
es.spineLength = sl(st <= tx(I) & et > tx(I) & et-st > 30);



es.animalTime = length(expt.gatherField('eti'))*dt;

es.numruns = length([expt.track.run]);
es.numturns = nnz(expt.gatherSubField('reorientation', 'numHS') > 0);
es.numpauses = nnz(expt.gatherSubField('reorientation', 'numHS') == 0);
