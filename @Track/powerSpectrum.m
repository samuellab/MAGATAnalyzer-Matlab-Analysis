function [ps, f] = powerSpectrum(track, quantityName, timeInterval, varargin)
% calculate the power spectrum of quantityName;
% function [ps, f] = powerSpectrum(track, quantityName, timeInterval, varargin)
%
% calculate the power spectrum of quantityName; uses a window timeInterval
%       long
% outputs:
%   PS: the power spectrum
%   F: the frequencies at which the power spectrum is defined
% inputs:
%   TRACK < Track
%   QUANTITYNAME: the name of the quantity - must be 1D; e.g.
%       'speed' is OK, 'vel' is not
%   TIMEINTERVAL: the size of the averaging window window
%   VARARGINL: optional parameter value pairs
%       'MeanSubtracted', true/false (if true, subtracts mean value before
%           computing psd)
%
%   see Spectrum; Spectrum/Welch; Spectrum/PSD

MeanSubtracted = true;
varargin = assignApplicable(varargin);

Hs = spectrum.welch('Hamming', timeInterval/track.dr.interpTime);
x = track.getDerivedQuantity(quantityName);
if (MeanSubtracted)
    x = x - repmat(mean(x,2),[1 length(x)]);
end
%size(x)
for j = 1:size(x,1)
    hpsd = Hs.psd(x(j,:), 'Fs', 1 / track.dr.interpTime, 'NormalizedFrequency', false);

    ps(j,:) = hpsd.Data;
    f(j,:) = hpsd.Frequencies/track.dr.interpTime;
end