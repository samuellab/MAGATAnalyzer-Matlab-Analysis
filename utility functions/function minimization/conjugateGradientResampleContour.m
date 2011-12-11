function x = conjugateGradient(fn, gradfn, x0, varargin) 
%function x = conjugateGradient(fn, gradfn, x0) 
%
%implements conjugate gradient descent, inspired by nr

maxiter = 200;
ftol = 1E-10;
resampleInterval = 10;
plotInterval = 0;
varargin = assignApplicable(varargin);

x = x0;
fp = fn(x);
xi = -gradfn(x);
g = xi;
h = g;
fv = zeros([1 maxiter]);
for its = 1:maxiter
    [x,fret] = linmin(fn, x, xi);
    fv(its) = fret;
    if (2*abs(fret - fp) < ftol*(abs(fret) + abs(fp) + eps))
        its
        figure(20); plot(1:its, fv(1:its));
        return;
    end
    fp = fret;
    oldx = x;
    if (mod(its, resampleInterval) == 0)
        x = resampleContour(x); %new, test
    end
    xi = gradfn(x); %new, test
    %{
    gg = sum(g.^2);
    dgg = sum((xi+g).*xi);
    if (gg == 0)
        disp('gg = 0');
        its
        return;
    end
    gam = dgg/gg;
    g = -xi;
    xi = g+gam*h;
    oldx = x;
    [x, g, xi] = resampleContour(x, g, xi);
    h = xi;
    %}
   %{
    dl = [0 sqrt(sum(diff(x(:,[1:end 1]),[],2)).^2)];
    npts = length(x);
    l = cumsum(dl);
    il = linspace(0, max(l)*(npts-1)/npts, npts);
    oldx = x;
    x = interp1(l, x(:,[1:end 1])', il)';
    g = interp1(l, g(:,[1:end 1])', il)';
    xi = interp1(l, xi(:,[1:end 1])', il)';
    h = xi;
    %}
    if (plotInterval > 0 && mod(its, plotInterval) == 0)
        figure(11);plot (oldx(1,:), oldx(2,:),'b.-', x(1,:), x(2,:), 'r.-'); hold on; quiver(x(1,:), x(2,:), xi(1,:), xi(2,:)); hold off; pause
    end
end
disp ('too many iters');


function [x,fval] = linmin(fn, x, dx)
myfn = @(a) fn(x + a*dx);
[a,b] = minbrak(myfn, 0, 0.01);
if (a > b)
    temp = a;
    a = b;
    b = temp;
end
%{
testvals = linspace(a,b,20);
for j = 1:length(testvals)
    fnvals(j) = myfn(testvals(j));
end
plot (testvals, fnvals, 'b.-'); pause
%}
[a,fval] = fminbnd(myfn, a, b);
x = x + a*dx;
