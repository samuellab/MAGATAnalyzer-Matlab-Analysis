function [meanPeriFreq, freqs, psd] = setDerivationRulesByPeristalsisFrequency (eset, varargin)
%function setDerivationRulesByPeristalsisFrequency (eset, varargin)

%homogeonize interpolation time
it = eset.gatherSubField('dr', 'interpTime');
interpTime = percentile(it, 0.05);
if (any (it ~= interpTime))
    warning ('ESET:ITIME', ['eset has non homogeneous interpolation times; updating to all have same - ' num2str(interpTime)]); 
    eset.evaluateTrackExpression(['track.dr.interpTime = ' num2str(interpTime, 10) ';']);
    eset.executeTrackFunction('recalculateDerivedQuantities'); 
end
if (interpTime >= 0.25)
    warning ('ESET:ITIME','Interpolation Time is too slow for likely peristalsis frequencies');
end

im = eset.gatherField('imid');
im = im(:,all(isfinite(im)));
sigma = min(1.5, 0.1/interpTime);
vm = sqrt(sum(deriv(im,sigma).^2))/interpTime;

%important: nuke large jumps between tracks!
vm(vm > percentile(vm, 0.98)) = percentile(vm, 0.98);

Hs = spectrum.welch('Hamming', 20/interpTime);
hpsd = Hs.psd(vm - mean(vm), 'Fs', 1 / interpTime, 'NormalizedFrequency', false); ps = hpsd.Data; f = hpsd.frequencies;
freqs = f;
psd = ps;

cutoffFrequency = 0.4;
ps = ps(f > cutoffFrequency);
f = f(f > cutoffFrequency);
[~,I] = max(ps);
periFreq = f(I);
if (nargout > 0)
    meanPeriFreq = periFreq;
end
%plot (f, ps, f(I), ps(I), 'r*');

if (I == 1 || I == length(ps))
    warning ('ESET:PERI','failed to detect a peristalsis peak, aborting');
    return;
end

smoothTime = 0.5/periFreq;
derivTime = smoothTime/2;

for j = 1:length(eset.expt)
    eset.expt(j).dr.interpTime = interpTime;
    eset.expt(j).dr.smoothTime = smoothTime;
    eset.expt(j).dr.derivTime = derivTime;
end
eset.evaluateTrackExpression(['track.dr.smoothTime = ' num2str(smoothTime, 10) ';']);
eset.evaluateTrackExpression(['track.dr.derivTime = ' num2str(derivTime, 10) ';']);
eset.evaluateTrackExpression(['track.so.minRunTime = ' num2str(2.5/periFreq, 10) ';']);
eset.evaluateTrackExpression(['track.so.smoothBodyTime = ' num2str(0.2/periFreq, 10) ';']);
eset.executeTrackFunction('recalculateDerivedQuantities');
eset.executeExperimentFunction('assignGlobalQuantities');