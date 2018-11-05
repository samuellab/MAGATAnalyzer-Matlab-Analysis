function [yv, yv_alpha] = histVarBW (xdata, xaxis, w, kernelType, causal, alpha)
%function yv = histVarBW (ydata, xaxis, w, kernelType, causal)
%
%adapted from ssvkernel by Hideaki Shimazaki 
% http://2000.jukuin.keio.ac.jp/shimazaki

existsAndDefault('kernelType', 'Gauss');
existsAndDefault('causal', false);
existsAndDefault('alpha', 0.95);
dx = median(diff(xaxis));

xdata = reshape(xdata,1,numel(xdata));

if (any(abs(diff(xaxis)-dx)/dx > 1e-6))
    warning ('expected evenly spaced xbins');
end
y_hist = histc(xdata, binEdgesFromCenters(xaxis)); y_hist = y_hist(1:end-1);
N = sum(y_hist);
L = length(xaxis);

yv = zeros(1,L);
yone = ones(size(y_hist));
nf = zeros(1,L);

if (ischar(kernelType))
    switch lower(kernelType)
        case 'gauss'
            kf = @Gauss;
        case 'laplace'
            kf = @Laplace;
        case 'cauchy'
            kf = @Cauchy;
        case 'boxcar'
            kf = @Boxcar;
    end
else
    kf = kernelType;
end
if causal
    fun = @(x,ww) kf(x,ww).*(x >= 0);
else
    fun = kf;
end
            
for k = 1:L
   yv(k) = sum( y_hist.*fun(xaxis(k)-xaxis,w(k)) );
   nf(k) = sum( yone.*fun(xaxis(k)-xaxis,w(k)) );
end
yv = yv ./ nf;
yv = yv / sum(yv) * N;

if (nargout < 2)
    return;
end
nbs = ceil(20/(1-alpha)); 
yb = zeros(nbs,length(xaxis));
for i = 1: nbs, %disp([i nbs])
    Nb = poissrnd(N); %choose a poisson distributed number of samples
    idx = ceil(rand(1,Nb)*N);
    xb = xdata(idx); %randomly select Nb samples from original data
    y_histb = histc(xb, binEdgesFromCenters(xaxis)); y_histb = y_histb(1:end-1); %make fine-scaled histogram
    
    
    ybuf = zeros(1,L);
    for k = 1:L
        ybuf(k) = sum( y_histb.*fun(xaxis(k)-xaxis,w(k)) );
        nf(k) = sum( yone.*fun(xaxis(k)-xaxis,w(k)) );
    end
    ybuf = ybuf ./ nf;
    ybuf = ybuf / sum(ybuf) * Nb;
    
    yb(i,:) = ybuf;
    
end

ybsort = sort(yb);
y95b = ybsort(floor((1-alpha)/2*nbs),:); %95 is legacy from this being a 95% confidence interval
y95u = ybsort(floor((1+alpha)/2*nbs),:);

yv_alpha = [y95b; y95u];



function y = Gauss(x,w) 
y = 1/sqrt(2*pi)./w.*exp(-x.^2/2./w.^2);

function y = Laplace(x,w)
y = 1./sqrt(2)./w.*exp(-sqrt(2)./w.*abs(x));

function y = Cauchy(x,w) 
y = 1./(pi*w.*(1+ (x./w).^2));

function y = Boxcar(x,w)
a = sqrt(12)*w;
%y = 1./a .* ( x < a/2 ) .* ( x > -a/2 );
%y = 1./a .* ( abs(x) < a/2 );
y = ones(size(x))./a; 
y(abs(x) > a/2) = 0; %speed optimization
