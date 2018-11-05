function [newctr, varargout] = resampleContour(ctr, varargin)
%resamples contour and any other associated quantities (varargins) to 
%be evenly spaced (associated quantities are resampled along the same
%indices as contour
%
closed = true;
varargin = assignApplicable(varargin);
if (closed)
    npts = length(ctr);
    dl = sqrt(sum(diff(ctr(:,[1 1:end 1]),[],2).^2));
    l = cumsum(dl);
    il = linspace(0, l(end)*(npts-1)/npts, npts);
    inds = interp1(cumsum(dl), 1:(npts + 1), il);

    %newctr = interp1(cumsum(dl), ctr(:, [1:end 1])', il,'linear')';
    newctr = interp1(ctr(:, [1:end 1])', inds, 'spline')';

    for j = 1:length(varargin)
        x = varargin{j};
        varargout{j} = interp1(x(:, [1:end 1])', inds, 'spline')';
    end
else
    npts = length(ctr);
    dl = sqrt(sum(diff(ctr(:,[1 1:end]),[],2).^2));
    l = cumsum(dl);
    il = linspace(0, l(end), npts);
    [l,I] = unique(l);
    inds = interp1(l, I, il);

    %newctr = interp1(cumsum(dl), ctr(:, [1:end 1])', il,'linear')';
    newctr = interp1(ctr', inds, 'spline')';

    for j = 1:length(varargin)
        x = varargin{j};
        varargout{j} = interp1(x', inds, 'spline')';
    end
end