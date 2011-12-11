function x = conjugateGradient(fn, gradfn, x0, varargin) 
%function x = conjugateGradient(fn, gradfn, x0) 
%
%implements conjugate gradient descent, inspired by nr
debug = false;
maxiter = 200;
ftol = 1E-6;
conjugate = true;
varargin = assignApplicable(varargin);

x = x0;
fp = fn(x);
xi = -gradfn(x);
g = xi;
h = g;
fv = zeros([1 maxiter]);
for its = 1:maxiter
    
    [x,fret] = linmin(fn, x, xi, false);
    fv(its) = fret;
    if (2*abs(fret - fp) < ftol*(abs(fret) + abs(fp) + eps))
        %its
        if (debug)
            figure(20); plot(1:its, fv(1:its));          
        end
        return;
    end
    fp = fret;
    xi = gradfn(x);
    if (false)
        figure(10); clf(10); plot (x(1,:), x(2,:), 'b.-'); hold on; quiver (x(1,:), x(2,:), -xi(1,:), -xi(2,:), 'g');
    end
    if (conjugate)     
        gg = sum(g.^2);
        dgg = sum((xi+g).*xi);
        if (gg == 0)
            disp('gg = 0');
            %its
            return;
        end
        gam = dgg/gg;
        g = -xi;
        xi = g+gam*h;
        h = xi;
         
    else
        xi = -xi;
    end
    if (false)
        quiver (x(1,:), x(2,:), xi(1,:), xi(2,:), 'r'); hold off;
        axis equal
        pause
    end
end
disp ('too many iters');
if (debug)
    figure(20); plot(1:its, fv(1:its));    
end

function [x,fval] = linmin(fn, x, dx, debug)
myfn = @(a) fn(x + a*dx);
[a,b] = minbrak(myfn, 0, 0.01);
if (a > b)
    temp = a;
    a = b;
    b = temp;
end
if (debug)
    figure(11); clf(11);
    aa = linspace(a, b, 20);
    for j = 1:length(aa)
        fv(j) = fn(x + aa(j)*dx);
    end
    plot (aa, fv);
end
%{
testvals = linspace(a,b,20);
for j = 1:length(testvals)
    fnvals(j) = myfn(testvals(j));
end
plot (testvals, fnvals, 'b.-'); pause
%}
[a,fval] = fminbnd(myfn, a, b);
if (debug)
    figure(10); clf(10); plot (x(1,:), x(2,:), 'b.-', x(1,:) - a*dx(1,:), x(2,:) - a*dx(2,:), 'r.-');
    pause;
end
x = x + a*dx;
