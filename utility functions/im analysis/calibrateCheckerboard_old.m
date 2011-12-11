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
%
%
%xc,yc are coordinates of corners
%structvar follows the rules of structvar, which are
%
%about structvar
%the first rule of structvar is you do not talk about structvar
%the second rule of structvar is you do not talk about structvar
%the third rule of structvar is that it is a Nx2 array of points,
%specifying rectangles
%the fourth rule of structvar is the lower left (in image coordinates; upper left in xy coords)
%corner comes first, then the other 3 points are specified in
%counterclockwise (in image coords; clockwise in xy coords) order
%the fifth rule of structvar is you do not talk about structvar
%
%structvar specifies the boundaries of light squares

xaxis = [];
yaxis = [];
flatten = false;
varargin = assignApplicable(varargin);
if (isempty(xaxis))
    xaxis = 1:size(im,2);
end
if (isempty(yaxis))
    yaxis = 1:size(im,1);
end

figure(10); clf(10); 
imagesc(xaxis, yaxis, double(im)); axis equal; colormap gray(256);   
title ('please select opposite corners of a single square in center');
   
while (1)
    [x,y] = getpts;
    if (length(x) < 2 || length(y) < 2) 
        figure(10); clf(10); 
        imagesc(xaxis, yaxis, double(im)); axis equal; colormap gray(256);   
        title ('please select two points');
        continue;
    end
    try 
        figure(10); clf(10); 
        imagesc(xaxis, yaxis, double(im)); axis equal; colormap gray(256);   
        title ('something messed up - try again');
        width = diff(x(1:2));
        height = diff(y(1:2));
    catch
        continue;
    end
    title ('thank you - doing spatial calibration, this may take a moment');
    pause(0.1);
    break;
end
if (flatten)
    checkerSize = ceil(max(width, height));
    intensityim = blurim(imdilate(blurim(im, 4), ones(ceil(checkerSize/2))), checkerSize/3);
    im = double(im) ./ double(intensityim);
    im(im > 2) = 2;
    figure(9); clf(9); imagesc(xaxis, yaxis, intensityim);
    figure(11); clf(11); imagesc(xaxis, yaxis, double(im)); axis equal; colormap jet;   
end
    
x = interp1(xaxis, 1:length(xaxis), x);
y = interp1(yaxis, 1:length(yaxis), y);

dx = abs(diff(x(1:2)));
dy = abs(diff(y(1:2)));



%look for corner locations by convolving with checker pattern

%construct kernel
ksize = [dy/1.5,  dx/1.5];
%make ksize odd
ksize = 2 * round(ksize/2)+1;

kernel = zeros(ksize);
kernel(1:floor(ksize(1)/2),1:floor(ksize(2)/2)) = 1;
kernel(ceil(ksize(1)/2 + 1):end,1:floor(ksize(2)/2)) = -1;
kernel(ceil(ksize(1)/2 + 1):end,ceil(ksize(2)/2 + 1):end) = 1;
kernel(1:floor(ksize(1)/2),ceil(ksize(2)/2 + 1):end) = -1;

kernel = kernel / sum(abs(kernel(:)));
%convolve & discard outer region 
cim = zeros(size(im));
cimall = conv2(double(im), kernel, 'same');
yinds = round(ksize(1)/3):round(size(cim,1) - ksize(1)/3);
xinds = round(ksize(2)/3):round(size(cim,2) - ksize(2)/3);
cim(yinds,xinds) = cimall(yinds,xinds);
%pcolor(xaxis,yaxis,cim); shading flat; colormap jet;

%look for local extrema

%do a 2nd deriv filter
sd = zeros(size(cim));
tempim = conv2(gaussKernel(ksize(1)/18), gaussKernel(ksize(2)/18), cim, 'same') - conv2(gaussKernel(ksize(1)/6), gaussKernel(ksize(2)/6), cim, 'same');
yinds = round(ksize(1)/2):round(size(cim,1) - ksize(1)/2);
xinds = round(ksize(2)/2):round(size(cim,2) - ksize(2)/2);
sd(yinds, xinds) = tempim(yinds,xinds);
sdlg = abs(sd) > percentile(abs(sd), 0.975);
figure(12); clf(12);
pcolor(xaxis, yaxis, sd); shading flat; colormap jet; axis equal
figure(13); clf(13);
pcolor(xaxis, yaxis, sd.*sdlg); shading flat; colormap jet; axis equal

%find the locations of the corners
stats = regionprops(imdilate(sdlg,strel('disk',10, 0)), abs(cim), 'WeightedCentroid');
loc = [stats.WeightedCentroid];
[xl,yl] = localmaxima(abs(cim), loc(1:2:end), loc(2:2:end),max(ksize)/4);
%return back to coordinate system specified by xaxis,yaxis
%the extra shift of 1/2 was determined empirically
xl = interp1(xaxis, xl + 1/2, 'linear');
yl = interp1(yaxis, yl + 1/2, 'linear');

xc = xl;
yc = yl;

%sort into rows
[yc,I] = sort(yc);
xc = xc(I);

jump = [0 find(diff(yc) > height * 0.5)];
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
width = median(dx);
height = median(diff(yavg));

x0 = min(xc);
y0 = min(yavg);
for j = 1:length(row)
    row(j).xreal = round((row(j).xc - x0)/width);
    row(j).yreal = repmat(round((yavg(j)-y0) / height),size(row(j).xreal));
end

camx = [row.xc];
camy = [row.yc];
realx = [row.xreal];
realy = [row.yreal];
figure(10);
hold on;
colormap gray;
plot (camx, camy, 'r.');
hold off
[realx, realy, camx, camy] = guessOutsideHull (realx, realy, camx, camy, -1:(size(im,2) + 1), -1:(size(im,1)+1));
realx = realx - floor(min(realx));
realy = realy - floor(min(realy));
camx = camx;
camy = camy;

title ('done');

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


    