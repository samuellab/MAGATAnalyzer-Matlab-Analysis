function [newctr, newctri] = marcsnake (alpha, beta, energyim, sigmas, oldctr)
%function newctr = marcsnake (alpha, beta, energyim, sigmas, oldctr)
%
%works to minimize the energy of the contour = alpha(dx^2) + beta ddx^2 +
%energyim(ctr)
%
%if sigmas is a vector, then we do successive minimizations with energyim
%blurred by sigmas(1), then sigmas(2),  etc.
%
%oldctr is the starting position.  the contour is assumed to be closed
%(derivatives are circular)
%
%ctr should be 2xNpts long

if (exist('sigmas', 'var') && length(sigmas) > 1)
    newctr = oldctr;
    for j = 1:length(sigmas)
        newctr = marcsnake(alpha, beta, energyim, sigmas(j), newctr);
       % pause
    end
    return;
end

existsAndDefault('sigmas', 1);
sigma = sigmas(1);

energyim = double(energyim);
[dxim, dyim] = imgradient(energyim, sigma);
energyim = blurim(energyim, sigma);

energyfn = @(x) snakeEnergy(x, alpha, beta, energyim, sigma);
gradfn = @(x) gradEnergy(x, alpha, beta, dxim, dyim, sigma);

[oldu, oldspr, oldbend, oldpot] = snakeEnergy(oldctr, alpha, beta, energyim, sigma);
[gx, gspr, gbend, gpot] = gradEnergy(oldctr, alpha, beta, dxim, dyim, sigma);
figure(3); clf();
plot(1:length(gx), gx(1,:), 1:length(gspr), gspr(1,:), 1:length(gbend), gbend(1,:), 1:length(gpot), gpot(1,:), 1:length(beta), beta);
max(sum(gbend.^2))
legend ('xd', 'spr', 'bend', 'pot', 'beta')

gnorm = sqrt(max(sum([gspr.^2 gbend.^2 gpot.^2], 2)))*5 / sqrt(sum(size(energyim).^2));
figure(1); clf();
pcolor(energyim); colormap gray; shading interp; hold on;
plot (oldctr(1,:), oldctr(2,:), 'r-','LineWidth',2); hold on;
quiver(oldctr(1,:), oldctr(2,:), -gx(1,:)./gnorm,-gx(2,:)./gnorm, 0,'r');
quiver(oldctr(1,:), oldctr(2,:), -gspr(1,:)./gnorm,-gspr(2,:)./gnorm, 0,'g');
quiver(oldctr(1,:), oldctr(2,:), -gbend(1,:)./gnorm,-gbend(2,:)./gnorm, 0,'c');
quiver(oldctr(1,:), oldctr(2,:), -gpot(1,:)./gnorm,-gpot(2,:)./gnorm,0,'y');
hold off

newctr = conjugateGradientResampleContour(energyfn, gradfn, oldctr, 'maxiter', 1000, 'ftol', 1E-12); 
%{
perturb = 0.5*[cos(linspace(0,5*pi,length(newctr)));sin(linspace(0,5*pi,length(newctr)))];
newctr = conjugateGradient(energyfn, gradfn, newctr+perturb, 'maxiter', 1000, 'ftol', 1E-16);
%}
figure(2); clf();
[newu, newspr, newbend, newpot] = snakeEnergy(newctr, alpha, beta, energyim, sigma);
[gx, gspr, gbend, gpot] = gradEnergy(newctr, alpha, beta, dxim, dyim, sigma);

%gnorm = sqrt(max(sum([gspr.^2 gbend.^2 gpot.^2], 2)))*5 / sqrt(sum(size(energyim).^2));
pcolor(energyim); colormap gray; shading interp; hold on;
plot (newctr(1,:), newctr(2,:), 'r.-','LineWidth',2,'MarkerSize',30); hold on;
quiver(newctr(1,:), newctr(2,:), -gx(1,:)./gnorm,-gx(2,:)./gnorm, 0,'r');
quiver(newctr(1,:), newctr(2,:), -gspr(1,:)./gnorm,-gspr(2,:)./gnorm, 0,'g');
quiver(newctr(1,:), newctr(2,:), -gbend(1,:)./gnorm,-gbend(2,:)./gnorm, 0,'c');
quiver(newctr(1,:), newctr(2,:), -gpot(1,:)./gnorm,-gpot(2,:)./gnorm,0,'y');
legend('im', 'ctr', 'grad', 'spr', 'bend', 'pot')
hold off
figure(4); clf();
plot(1:length(gx), gx(1,:), 1:length(gspr), gspr(1,:), 1:length(gbend), gbend(1,:), 1:length(gpot), gpot(1,:));
figure(5);
plot(newctr(1,:), newctr(2,:), newctr(1,:)-10*gx(1,:), newctr(2,:)-10*gx(2,:));
energyfn(newctr)
dval = linspace(-1, 1, 100);
for j = 1:length(dval)
    [u, uspr, ubend, upot] = snakeEnergy(newctr+dval(j)*gx, alpha, beta, energyim, sigma);
    uu(j) = u;
    uuspr(j) = uspr;
    uubend(j) = ubend;
    uupot(j) = upot;
    
end
figure(6);
    plot (dval, uu-newu, dval, uuspr-newspr, dval, uubend-newbend, dval, uupot-newpot); 
    legend ('u', 'uspr', 'ubend', 'upot');
newu - oldu

function [u,uspr,ubend,upot] = snakeEnergy(ctr, alpha, beta, im, sigma)
upot = potenergy(ctr, im);
uspr = springenergy(ctr)*alpha;
ubend =bendenergy(ctr)*beta;
u = uspr + ubend + upot;
return;

n = length(ctr);
%dx = diff(ctr(:, [end 1:end]),[],2);
dx = deriv(ctr, sigma, 'padtype', 'circular');
dl = sqrt(sum(dx.^2));
L = sum(dl);
uspr = L*alpha*sum(1./dl)/(n^2);
that = dx./[dl;dl];
%cs = sum((diff(that(:,[end 1:end]), [],2)).^2)./(dl.^2); %curvature squared
cs = sum(deriv(that, sigma, 'padtype', 'circular').^2)./(dl.^2);
ubend = sum(beta.*cs)*L/n;
upot = sum((dl/L).*interp2(im, ctr(1,:), ctr(2,:), '*linear', 100*max(abs(im(:)))));
u = uspr + ubend + upot;

function [dx, gspr, gbend, gpot] = gradEnergy(ctr, alpha, beta,  xderivim, yderivim, sigma)
gpot = gradpotenergy(ctr, xderivim, yderivim);
gspr = gradspringenergy(ctr)*alpha;
gbend = gradbendenergy(ctr)*beta;

dx = gspr+gbend +gpot;
magdx = sqrt(sum(sum(dx.^2))) + eps;
gpot = gpot/magdx;
gspr = gspr/magdx;
gbend = gbend/magdx;
dx = dx/magdx;

return;

n = length(ctr);
%dx = diff(ctr(:, [end 1:end]),[],2);
dx = deriv(ctr, sigma, 'padtype', 'circular');
dl = sqrt(sum(dx.^2));
L = sum(dl);
dll = [dl;dl];
%{
L = sum(dl);
dL = sum(dl)/n;
dll = [dl;dl];
ddx = diff(ctr(:, [end 1:end 1]), 2, 2);
d4x = diff(ctr(:, [end-1:end 1:end 1:2]), 4, 2);
%}
d4x = dx;
for j = 1:3
    d4x = deriv(d4x, sigma, 'padtype', 'circular');
end
dx1 = dx./([dl.^3;dl.^3]);
gspr = alpha * deriv(dx1, sigma, 'padtype', 'circular');
%{
dx1 = -dx./([dl.^3;dl.^3]);
dx2 = -dx1(:,[2:end 1]);
gspr = alpha * (dx1 + dx2);
%}
gbend = -L*[beta.*d4x(1,:);beta.*d4x(2,:)]./(dll.^4);%./(dll.^4);

fx = interp2(xderivim, ctr(1,:), ctr(2,:), '*linear', NaN);
fy = interp2(yderivim, ctr(1,:), ctr(2,:), '*linear', NaN);

%nnz(~isfinite(fx))
fx(~isfinite(fx)) = mean(ctr(1,:)) - fx(~isfinite(fx));
fy(~isfinite(fy)) = mean(ctr(2,:)) - fy(~isfinite(fy));
gpot = [fx;fy];
dx = gspr+gbend +gpot;