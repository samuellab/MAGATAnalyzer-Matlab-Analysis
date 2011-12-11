function [realx, realy, camx, camy] = calibrateCheckerboard(im, varargin)
%function [realx, realy, camx, camy] = calibrateCheckerboard(im, varargin)
%varargin)
%   finds the corners of a checkerboard image, with some help
%   xaxis, yaxis are the axes for the image (i.e. pcolor(axis,yaxis,im)
%   makes sense)
%   im must be a grayscale image (size nrowsxncolsx1)
%   checkers must be roughly aligned to axes
%optional key/value paris
%"xaxis" x-axis of im
%"yaxis" y-axis of im
%"flatten", true/false -- whether to flatten out irregularities in
%illumination
%"xinds", xi "yinds",yi only consider these indices when calibrating checkerboard
%"flipy", [false]/true, if true, make low camera y values high real values 
%"flipx", [false]/true, if true, make low camera x values high real values
%realx,realy camx,camy mark the correspondence between real points (1
%square = 1 unit)
%"affine", [false]/true, if true, don't assume that checkerboard is almost
%regular when calibrating

xaxis = [];
yaxis = [];
xinds = [];
yinds = [];
flipy = false;
flipx = false;
flatten = false;
harrissigma = 3;
affine = false;
varargin = assignApplicable(varargin);
origimsize = size(im);
if (~isempty(xinds))
    xaxis = xinds;
    if (isempty(yinds))
        im = im(:,xinds);       
    else
        im = im(yinds, xinds);
        yaxis = yinds;
    end
else
    if (~isempty(yinds))
        im = im(yinds, :);
        yaxis = yinds;
    end
end

if (isempty(xaxis))
    xaxis = 1:size(im,2);
end
if (isempty(yaxis))
    yaxis = 1:size(im,1);
end

if (flatten)
    [~, cpts] = myharris(im, harrissigma);
    checkerSize = sqrt(size(im,1)*size(im,2)/length(cpts));
    intensityim = blurim(imdilate(blurim(im, 4), ones(ceil(checkerSize/1.5))), checkerSize/2);
    im = double(im) ./ double(intensityim);
    im(im > 2) = 2;
%     figure(9); clf(9); imagesc(xaxis, yaxis, intensityim);
%     figure(11); clf(11); imagesc(xaxis, yaxis, double(im)); axis equal; colormap jet;   
end

[cim, cpts] = myharris(im, harrissigma);

%whos cpts

% figure(21); clf(21);
% imagesc(xaxis, yaxis, cim); axis equal; colormap jet

xc = interp1(xaxis, cpts(1,:));
yc = interp1(yaxis, cpts(2,:));

%sort into rows

[yc,I] = sort(yc);
xc = xc(I);

[~,c] = kmeans(diff(yc), 2, 'start', [0; max(diff(yc))]);
height = max(c);


jump = unique([0 find(diff(yc) > height * 0.5) length(yc)]);
for j = 1:(length(jump) - 1)
    xx = xc((jump(j) + 1):(jump(j+1)));
    yy = yc((jump(j) + 1):(jump(j+1)));
    [xx,I] = sort(xx);
    yy = yy(I);
    row(j).xc = xx; %#ok<AGROW>
    row(j).yc = yy; %#ok<AGROW>
    dx(j) = median(diff(row(j).xc));
    yavg(j) = mean(row(j).yc);
end
width = median(dx(isfinite(dx)));
height = median(diff(yavg));

x0 = median(mod(xc,width));
y0 = median(mod(yc,height));

xr = (xc - x0)/ width;
yr = (yc - y0)/ height;

if (affine)
    [xc, yc, xr, yr] = affineFit(xc, yc, xr, yr);
    [~,I] = sort(xc);
    [~,J] = sort(yr(I));
    w = median(diff(xc(I(J))))
    [~,I] = sort(yc);
    [~,J] = sort(xr(I));
    h = median(diff(yc(I(J))))
    
else
        %remove any points that aren't within 10% of the checker size of where the
    %corner would be if the grid were perfect
    valid = (abs(xr - round(xr)) < width/10) & (abs(yr - round(yr)) < height/10);
    xrdub = xr(valid);
    yrdub = yr(valid);
    xc = xc(valid);
    yc = yc(valid);


    xr = round(xrdub);
    yr = round(yrdub);

    %if two or more camera points collapse to the same real point, choose the
    %closest to the actual real point
    valid = true(size(xr));
    for j = 1:length(xr)
        if (~valid(j))
            continue
        end
        same = find((xr((j+1):end) == xr(j)) & (yr((j+1):end) == yr(j)));
        if (isempty(same))
            continue
        end
        inds = [j same+j];
        [~,I] = min((xrdub(inds)-xr(j)).^2 + (yrdub(inds)-yr(j)).^2);
        valid(inds) = false;
        valid(inds(I)) = true;
    end

    xr = xr(valid);
    yr = yr(valid);
    xc = xc(valid);
    yc = yc(valid);

    %remove points within w,h of the edge
    h = ceil(height/5);
    w = ceil(width/5);
    sz = origimsize;
    valid = floor(xc) > w & ceil(xc) < sz(2) - w & floor(yc) > h & ceil(yc) < sz(1) - h;
    xr = xr(valid);
    yr = yr(valid);
    xc = xc(valid);
    yc = yc(valid);

    %remove points that aren't at the intersection of two checkers
    %{
    [idx,c] = kmeans(double(im(:)), 2, 'start', double([min(im(:)); max(im(:))]));
    if (c(2) > c(1))
        high = 2;
    else
        high = 1;
    end
    bwim = reshape(idx == high, size(im));
    %}
end

valid = true(size(xc));
for j = 1:length(xc)
        
    xim = interp1(xaxis, 1:length(xaxis), [xc(j)-w, xc(j)+w], 'nearest','extrap');
    yim = interp1(yaxis, 1:length(yaxis), [yc(j)-h, yc(j)+h], 'nearest','extrap');
    xi = min(xim):max(xim);
    yi = min(yim):max(yim);
%    xi = round(xc(j)-w):round(xc(j)+w);
%    yi = round(yc(j)-h):round(yc(j)+h);
    imdat = double(im(yi,xi));
    idx = kmeans(imdat(:), 2,'start',[min(imdat(:)); max(imdat(:))]);
 %   idxmean(j) = mean(idx);
    valid(j) = abs(mean(idx)-1.5) < 0.1;
end
%figure(1); clf; plot(idxmean);
%figure(2); clf; hist(idxmean, 1:0.01:2);
%figure(3); plotColorLine(xc, yc, idxmean);

xr = xr(valid);
yr = yr(valid);
xc = xc(valid);
yc = yc(valid);



[xr,I] = sort(xr, 'ascend');
yr = yr(I);
xc = xc(I);
yc = yc(I);

[yr,I] = sort(yr, 'ascend');
xr = xr(I);
xc = xc(I);
yc = yc(I);


    

camx = xc;
camy = yc;
realx = xr;
realy = yr;
if (flipy)
    realy = -realy;
end
if (flipx)
    realx = -realx;
end
%{
for j = 1:length(row)
    row(j).xreal = round((row(j).xc - x0)/width);
    row(j).yreal = repmat(round((yavg(j)-y0) / height),size(row(j).xreal));
end

camx = [row.xc];
camy = [row.yc];
realx = [row.xreal];
realy = [row.yreal];
%}
figure(10); clf(10); 
imagesc(xaxis, yaxis, double(im)); axis equal; colormap gray(256);   
hold on
plot (camx, camy, 'r.');
hold off
[realx, realy, camx, camy] = guessOutsideHull (realx, realy, camx, camy, -1:(origimsize(2) + 1), -1:(origimsize(1)+1));
realx = realx - floor(min(realx));
realy = realy - floor(min(realy));
%camx = camx;
%camy = camy;

title ('result of checkerboard calibration -- red dots should lie on corners');

function [x,y] = localmaxima(im, x, y, nsize, maxiters, tol)
%function [x,y] = localmaxima(im, x, y, maxiters, tol)
%uses steepest ascent to move to local maxima;  
%this is really only a good idea if you're already close

debug = false;

existsAndDefault('nsize', 10);
existsAndDefault('maxiters', 100);
existsAndDefault('tol', 0.001);

xd = conv2(gaussKernel(nsize), dgausskernel(nsize), im, 'same');
yd = conv2(dgausskernel(nsize), gaussKernel(nsize), im, 'same');
ds = sqrt(xd.^2 + yd.^2);
xd = xd./ds;
yd = yd./ds;

xx = 1:size(im,2);
yy = 1:size(im,1);

xnew = x;
ynew = y;
nv = interp2(xx,yy,im,xnew,ynew);
if (debug)
    figure(1); clf(1); imagesc(im); 
    hold on
end
for j = 1:maxiters
    ov = nv;
    xdd = interp2(xx,yy,xd,x,y,'*linear');
    ydd = interp2(xx,yy,yd,x,y,'*linear');
    for k = 1:length(x)
        xline = x(k) + xdd(k)*(0:0.01:1)*nsize;
        yline = y(k) + ydd(k)*(0:0.01:1)*nsize;
        valsonline = interp2(xx,yy,im,xline,yline,'*linear');
        [blah,I] = max(valsonline);
        xnew(k) = xline(I);
        ynew(k) = yline(I);
    end
    if (debug)
        for k = 1:length(x)
            plot ([x(k) xnew(k)], [y(k) ynew(k)], 'b-', 'LineWidth',2);
        end
    end
    x = xnew;
    y = ynew;
    nv = interp2(xx,yy,im,xnew,ynew);
    if all((nv - ov)./ov < tol)
        return
    end
end


    