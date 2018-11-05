function newctr = marcsnake (alpha, beta, imxaxis, imyaxis, energyim, sigmas, oldctr)
%function newctr = marcsnake (alpha, beta, imxaxis, imyaxis, energyim, sigmas, oldctr)
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

existsAndDefault('sigmas', 1);
existsAndDefault('imxaxis', 1:size(energyim, 2));
existsAndDefault('imyaxis', 1:size(energyim, 1));
if (length(sigmas) > 1)
    newctr = oldctr;
    for j = 1:length(sigmas)
        newctr = marcsnake(alpha, beta, imxaxis, imyaxis, energyim, sigmas, oldctr);
    end
    return;
end
sigma = sigmas(1);

energyim = double(energyim);
[dxim, dyim] = imgradient(energyim, sigma);
energyim = blurim(energyim, sigma);
dx = diff(imxaxis(1:2)); dy = diff(imyaxis(1:2));
dxim = dxim./dx;
dyim = dyim./dy;
energyfn = @(x) snakeEnergy(x, alpha, beta, imxaxis, imyaxis, energyim);
gradfn = @(x) gradEnergy(x, alpha, beta, imxaxis, imyaxis, dxim, dyim);

[oldu, oldspr, oldbend, oldpot] = snakeEnergy(oldctr, alpha, beta, imxaxis, imyaxis, energyim)
[gx, gspr, gbend, gpot] = gradEnergy(oldctr, alpha, beta, imxaxis, imyaxis, dxim, dyim);
plot(1:length(gx), gx(1,:), 1:length(gspr), gspr(1,:), 1:length(gbend), gbend(1,:), 1:length(gpot), gpot(1,:));
max(sum(gbend.^2))
legend ('xd', 'spr', 'bend', 'pot')
return
gnorm = sqrt(max(sum([gspr.^2 gbend.^2 gpot.^2], 2)))*5 / sqrt(diff(imxaxis([1 end])).^2 + diff(imyaxis([1 end])).^2);
figure(1); clf();
pcolor(imxaxis, imyaxis, energyim); colormap gray; shading interp; hold on;
plot (oldctr(1,:), oldctr(2,:), 'r-','LineWidth',2); hold on;
quiver(oldctr(1,:), oldctr(2,:), -gx(1,:)./gnorm,-gx(2,:)./gnorm, 0,'r');
quiver(oldctr(1,:), oldctr(2,:), -gspr(1,:)./gnorm,-gspr(2,:)./gnorm, 0,'g');
quiver(oldctr(1,:), oldctr(2,:), -gbend(1,:)./gnorm,-gbend(2,:)./gnorm, 0,'c');
quiver(oldctr(1,:), oldctr(2,:), -gpot(1,:)./gnorm,-gpot(2,:)./gnorm,0,'y');
hold off

newctr = conjugateGradient(energyfn, gradfn, oldctr, 'maxiter', 1000);
figure(2); clf();
[newu, newspr, newbend, newpot] = snakeEnergy(newctr, alpha, beta, imxaxis, imyaxis, energyim)
[gx, gspr, gbend, gpot] = gradEnergy(newctr, alpha, beta, imxaxis, imyaxis, dxim, dyim);

gnorm = sqrt(max(sum([gspr.^2 gbend.^2 gpot.^2], 2)))*5 / sqrt(diff(imxaxis([1 end])).^2 + diff(imyaxis([1 end])).^2);
pcolor(imxaxis, imyaxis, energyim); colormap gray; shading interp; hold on;
plot (newctr(1,:), newctr(2,:), 'r-','LineWidth',2); hold on;
quiver(newctr(1,:), newctr(2,:), -gx(1,:)./gnorm,-gx(2,:)./gnorm, 0,'r');
quiver(newctr(1,:), newctr(2,:), -gspr(1,:)./gnorm,-gspr(2,:)./gnorm, 0,'g');
quiver(newctr(1,:), newctr(2,:), -gbend(1,:)./gnorm,-gbend(2,:)./gnorm, 0,'c');
quiver(newctr(1,:), newctr(2,:), -gpot(1,:)./gnorm,-gpot(2,:)./gnorm,0,'y');
hold off
newu - oldu

function [u,uspr,ubend,upot] = snakeEnergy(ctr, alpha, beta, imxaxis, imyaxis, im)

n = length(ctr);
dx = diff(ctr(:, [end 1:end]),[],2);
dl = sqrt(sum(dx.^2));
L = sum(dl);

ddx  = diff(ctr(:, [end-1:end 1:end]), 2, 2);

if (false)
    uspr = 1/n*alpha/2*sum(sum(dx.^2)/dL^2);
    ubend = 1/n*beta/2*sum(sum(ddx.^2)/dL^4);
    upot = sum(1/n*interp2(imxaxis, imyaxis, im, ctr(1,:), ctr(2,:), '*linear', 100*max(abs(im(:)))));
else
    uspr = L*alpha/2;
    %uspr = 1/n*alpha/2*sum(sum(dx.^2)/dL^2);
    that = dx./[dl;dl];
    cs = sum((diff(that(:,[end 1:end]), [],2)).^2)./(dl.^2); %curvature squared
    ubend = beta*sum(cs.*dl);
    %ubend = beta/2*sum(sum(ddx.^2)./(dl.^3));
    upot = sum(dl.*interp2(imxaxis, imyaxis, im, ctr(1,:), ctr(2,:), '*linear', 100*max(abs(im(:)))));
end
u = uspr + ubend + upot;

function [dx, gspr, gbend, gpot] = gradEnergy(ctr, alpha, beta, imxaxis, imyaxis, xderivim, yderivim)
n = length(ctr);
dx = diff(ctr(:, [end 1:end]),[],2);
dl = sqrt(sum(dx.^2));
dL = sum(dl)/n;
dll = [dl;dl];
ddx = diff(ctr(:, [end 1:end 1]), 2, 2);
d4x = diff(ctr(:, [end-1:end 1:end 1:2]), 4, 2);
if (true)
    gspr = -alpha*ddx/(dL.^2);%./(dll.^2);
    gbend = beta*d4x/(dL.^4);%./(dll.^4);
    fx = interp2(imxaxis, imyaxis, xderivim, ctr(1,:), ctr(2,:), '*linear', NaN);
    fy = interp2(imxaxis, imyaxis, yderivim, ctr(1,:), ctr(2,:), '*linear', NaN);
else
    gspr = -alpha*ddx./(dll.^2);
    gbend = beta*d4x./(dll.^4);
    fx = interp2(imxaxis, imyaxis, xderivim, ctr(1,:), ctr(2,:), '*linear', NaN);
    fy = interp2(imxaxis, imyaxis, yderivim, ctr(1,:), ctr(2,:), '*linear', NaN);
end
%nnz(~isfinite(fx))
fx(~isfinite(fx)) = mean(ctr(1,:)) - fx(~isfinite(fx));
fy(~isfinite(fy)) = mean(ctr(2,:)) - fy(~isfinite(fy));
gpot = [fx;fy];
dx = gspr+gbend +gpot;