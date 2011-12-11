function [ps, f] = powerSpectrum(expt, quantityName, timeInterval, varargin)
%produces the power spectrum for a qiven measured/derived quantity
%function [ps, f] = powerSpectrum(expt, quantityName, timeInterval,
%varargin)
%
%outputs: 
%PS: the power spectrum
%F: the frequencies at which the power spectrum is determined
%inputs:
%EXPT: a member of the experiment class
%TIME INTERVAL: the size of the window for averaging purposes
%   a smaller window has greater averaging but less frequency resolution
%VARARGIN: any parameters that can be passed to TRACK/POWERSPECTRUM
w = [expt.track.npts];
w = w/sum(w);
[ps,f] = expt.track(1).powerSpectrum(quantityName, timeInterval, varargin);

ps = zeros(size(ps));

for j = 1:length(expt.track)
    ps = ps + w(j)*expt.track(j).powerSpectrum(quantityName, timeInterval, varargin);
end


    