function x = chooseFromPDF (xaxis, pdf, n)
%function x = chooseFromPDF (xaxis, pdf, n)
%given a probability density function p = pdf(x) defined on xaxis, chooses
%an x within xaxis s.t. x is distributed according to pdf(x)
%returns n values
if (~exist('n', 'var') || isempty(n))
    n = 1;
end

pdf = pdf ./ sum(pdf);
%pdf(pdf <= 0) = 1E-6; %eliminate zero (should not have negative but get rid of those two) probabilities, so that cdf is strictly increasing

%pdf = pdf ./ sum(pdf);
cdf = cumsum(pdf);
[cdf,I] = unique(cdf);

x = interp1 (cdf, xaxis(I), rand([1 n]),'linear','extrap');



