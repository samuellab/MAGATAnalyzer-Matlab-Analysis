function [model, ctr] = fitBezierToImage (im, model0)
%c = fitLegendrePolynomialToImage (im, n)
im = double(im);

op = optimset('fminunc');
op.Display = 'off';
op.LargeScale = 'off';
doplots = true;

scale = 1;
im2 = imresize(im,scale);

op.DiffMinChange = 0.5;
model0.spineCpts = model0.spineCpts*scale;
model0.bodyWidth = model0.bodyWidth*scale;
myfun = @(q) energyfun(q, blurim(im2,2),model0);
if (doplots)
    op.PlotFcns = @(x, optimValues, state) myPlotFun(x,blurim(im2,2),model0);
end
x0 = modelTo1D(model0);
x0 = x0(1:8);
x1 = fminunc(myfun, x0, op);
model = model0;
model.spineCpts = [x1(1:4); x1(5:8)];
% 
% x0 = [model.spineCpts(1,:) model.spineCpts(2,:)];
% myfun = @(q) energyfun(q, im2,model);
% op.PlotFcns = @(x, optimValues, state) myPlotFun(x,im2,model);
% 
% op.DiffMinChange = 0.5;
% x1 = fminunc(myfun, x0, op);
% model.spineCpts = [x1(1:4); x1(5:8)];
op.DiffMinChange = 0.1;
myfun = @(q) energyfun(q, blurim(im2,2));
if (doplots)
    op.PlotFcns = @(x, optimValues, state) myPlotFun(x,blurim(im2,2));
end
x0 = modelTo1D(model);
x1 = fminunc(myfun, x0, op);
model = modelFrom1D(x1);


op.DiffMinChange = 0.5;

myfun = @(q) energyfun(q, im2,model);
if (doplots)
    op.PlotFcns = @(x, optimValues, state) myPlotFun(x,im2,model);
end
x0 = modelTo1D(model);
x0 = x0(1:8);
x1 = fminunc(myfun, x0, op);
model.spineCpts = [x1(1:4); x1(5:8)];
% 
% x0 = [model.spineCpts(1,:) model.spineCpts(2,:)];
% myfun = @(q) energyfun(q, im2,model);
% op.PlotFcns = @(x, optimValues, state) myPlotFun(x,im2,model);
% 
% op.DiffMinChange = 0.5;
% x1 = fminunc(myfun, x0, op);
% model.spineCpts = [x1(1:4); x1(5:8)];
op.DiffMinChange = 0.25;
myfun = @(q) energyfun(q, im2);
if (doplots)
    op.PlotFcns = @(x, optimValues, state) myPlotFun(x,im2);
end
x0 = modelTo1D(model);
x1 = fminunc(myfun, x0, op);
model = modelFrom1D(x1);



model.spineCpts = model.spineCpts/scale;
model.bodyWidth = model.bodyWidth/scale;
ctr = maggotContourFromControlPoints(model);
eic = energyInContour(im, ctr(1,:), ctr(2,:))
efie = energyFromImageEdge(ctr,im)
ecl = energyFromContourLength(x1)
ers = energyFromRigidSpine(x1)
efss = energyFromSpringSpine(x1)

end

function x = modelTo1D(model)
    x = [model.spineCpts(1,:) model.spineCpts(2,:) model.bodyWidth model.bodyTaper model.headalpha model.tailalpha];
end

function model = updateSpinePoints(x,model)
    model.spineCpts = [x(1:4); x(5:8)];
end
% function spinePointEnergyFun(q,im,model)
%     model = updateSpinePoints(q,model);
%     ctr = maggotContourFromControlPoints(model);
%     en = energyInContour(im, ctr(1,:), ctr(2,:))...
%         + 100*sum(sqrt(sum(diff(ctr,[],2).^2))) ...
%         + 100*sum(sum(diff(model.spineCpts,[],2).^2));
% end    

function model = modelFrom1D(x)
    model.spineCpts = [x(1:4); x(5:8)];
    model.bodyWidth = x(9);
    model.bodyTaper = x(10);
    model.headalpha = x(11);
    model.tailalpha = x(12);
end

function en = energyfun (q, im, model)
    if (nargin < 3)
        model = modelFrom1D(q);
    else
        model = updateSpinePoints(q,model);
    end
    ctr = maggotContourFromControlPoints(model);
    en = energyInContour(im, ctr(1,:), ctr(2,:))...
         + sum(sqrt(sum(diff(ctr,[],2).^2)))...
     + energyFromRigidSpine(q,model) ...
     + energyFromSpringSpine(q,model) + ...
     1E5*(1-model.tailalpha)^4 + 1E5*(2-model.headalpha)^4 + ...
     1E6*(0.8-model.bodyTaper)^4;
%         + sum(sum(diff(model.spineCpts,[],2).^2));
%         + energyFromImageEdge(ctr,im) ...

end    
function en = energyFromImageEdge(ctr, im)
%     if (nargin < 3)
%         model = modelFrom1D(q);
%     else
%         model = updateSpinePoints(q,model);
%     end
    edgim = conv2(dgausskernel(1), gaussKernel(1), im).^2 + conv2(gaussKernel(1), dgausskernel(1), im).^2;
    en = -sum(interp2(edgim,ctr(1,:),ctr(2,:),'*linear'))/sum(sqrt(sum(diff(ctr,[],2).^2)));
end

function en = energyFromContourLength(q,model)
    if (nargin < 2)
        model = modelFrom1D(q);
    else
        model = updateSpinePoints(q,model);
    end
     ctr = maggotContourFromControlPoints(model);
    en = sum(sqrt(sum(diff(ctr,[],2).^2)));
end
% function en = energyFromRigidSpine(q,model)
%     if (nargin < 2)
%         model = modelFrom1D(q);
%     else
%         model = updateSpinePoints(q,model);
%     end
%     sp = bezierFromControlPoints(model.spineCpts, 100);
%     dx = deriv(sp,1);
%     ddx = deriv(dx,1);
%     ds = sqrt(sum(dx.^2));
% %     rc = ds.^3./(dx(1,:).*ddx(2,:) - dx(2,:).*ddx(1,:));
% %     en = sum(ds./rc);
%    en = sum(abs((dx(1,:).*ddx(2,:) - dx(2,:).*ddx(1,:)))./ds.^2);
% end

function en = energyFromRigidSpine(q,model)
    if (nargin < 2)
        model = modelFrom1D(q);
    else
        model = updateSpinePoints(q,model);
    end
    t = diff(model.spineCpts,[], 2);
    s = sqrt(sum(t.^2));
    t = t./([1;1]*s);
    en = 1E5*sum(abs((t(1,1:(end-1)).*t(2,2:end) - t(2,1:(end-1)).*t(1,2:end))));
    
end
function en = energyFromSpringSpine(q,model)
    if (nargin < 2)
        model = modelFrom1D(q);
    else
        model = updateSpinePoints(q,model);
    end
    t = diff(model.spineCpts,[], 2);
    en = sum(sum(t.^2).^2) + 100*sum(sqrt(sum(t.^2)))*sum(1./(sum(t.^2)));
    
end

function stop = myPlotFun (q, im, model)
    imagesc(im);axis equal; hold on
    if (nargin < 3)
        model = modelFrom1D(q);
    else
        model = updateSpinePoints(q,model);
    end
    ctr = maggotContourFromControlPoints(model);
    plot(ctr(1,:), ctr(2,:), 'r-', model.spineCpts(1,:), model.spineCpts(2,:), 'k*--'); 
    stop = false;
end
function en = energyInContour (im, ctrx,ctry)
 %   bwim = traceBoundarySubPixel([ctrx;ctry],[1 1], size(im));

  %  en = sum((im(:) - mean(im(:).*bwim(:))).^2.*bwim(:)) ...
  %      + sum((im(:) - mean(im(:).*(1-bwim(:)))).^2.*(1-bwim(:)));
%    imagesc(bwim); axis equal; colorbar vert; pause(0.01);

bwim = poly2mask(ctrx, ctry, size(im,1), size(im,2));
%imagesc(im); axis equal; colorbar vert; hold on; plot(ctrx, ctry, 'r-'); hold off; pause(0.001);
en = sum((im(bwim)-mean(im(bwim))).^2) + sum((im(~bwim)-mean(im(~bwim))).^2); 

end