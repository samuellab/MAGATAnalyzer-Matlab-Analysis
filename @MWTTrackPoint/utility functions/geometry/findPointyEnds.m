function [inds,cv] = findPointyEnds (cpts, varargin)
%function inds = findPointyEnds (cpts)

sigma = length(cpts)/10.0;
debug = false;
varargin = assignApplicable(varargin);
ncp = length(cpts);

%{
cptspad = [cpts cpts];

cpl = lowpass1D(cptspad,sigma);
cpl = circshift(cpl(:,floor(ncp/2) + (1:ncp)), [0 -floor(ncp/2)]);

dx = 0.5*(diff(cpl(:, [end 1:end]),[],2) + diff(cpl(:, [1:end 1]),[],2));
dl = sqrt(sum(dx.^2));
that = dx./[dl;dl];
cv = sum((diff(that(:,[end 1:end]), [],2)))./(dl); %curvature
%}
cpl = lowpass1D(cpts,sigma,'padType', 'circular');
v = deriv(cpl, 1, 'padType', 'circular');
a = deriv(v, 1, 'padType', 'circular');

cv = (v(1,:).*a(2,:) - a(1,:).*v(2,:))./(sum(v.^2).^(1.5));

cv = cv .* sign(sum(cv));
[cvmax, maxI] = max(cv);

inds = 1:length(cv);
inddist = min(mod(inds-maxI+ncp,ncp), (mod (maxI - inds +ncp,ncp)));
%cvcorr = cv./cvmax + inddist/ncp;
cvcorr = cv;
[~,I] = sort(cvcorr,'descend');

ind1 = maxI;
ind2 = find(inddist(I) > ncp/4, 1, 'first');

inds = [ind1 I(ind2)];
if (debug)
    figure(1);
    
%cpl = circshift(cpl(:,floor(ncp/2) + (1:ncp)), [0 -floor(ncp/2)]);
    plot (cpl(1,:), cpl(2,:), 'k-'); hold on
    plot ([cpl(1,:); cpts(1,:)], [cpl(2,:); cpts(2,:)], 'g-', cpts(1,:), cpts(2,:), 'b-');
    plot (cpts(1,inds), cpts(2,inds), 'ro');
    plotColorLine(cpts(1,:), cpts(2,:), cv, jet); hold off
    axis equal
    figure(2);
    plot (cvcorr)
end
