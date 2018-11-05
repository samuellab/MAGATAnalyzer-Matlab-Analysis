function drawTrackImage (pt, camcalinfo, varargin)
%function drawTrackImage (pt, camcalinfo, varargin)
%@ImTrackPoint
%
%optional arguments, with defaults
%fid = [];
%scale = 1;
%cropSize = [];
%Axes = [];
%patchOptions = {};
%underlayImData = []; underlayImData.x, underlayImData.y, underlayImData.im
%showCov = false;
%
%clipToContour = false; if true, restrict image to region in and near the
%                      contour; only active if pt has a contour or 
%                      if contour is passed as an extra argument
%contour = (pt.contour) or [] ; a 2xN contour of the animal

fid = [];
existsAndDefault('camcalinfo', []);
scale = 1;
cropSize = [];
Axes = [];
underlayImData = [];
patchOptions = {};
showCov = false;
clipToContour = false;
underlayScale = -1;
if (any(strcmp('contour', properties(pt))))
    contour = double(pt.contour);
else
    contour = [];
end
varargin = assignApplicable(varargin);

if (isempty(Axes))
    Axes = gca;
end

if (isempty(pt.imData))
    if (~isempty(fid) && fid > 0)
        try
            fseek(fid, pt.locInFile, -1);
            pt2 = pt.fromFile(fid, true, true, camcalinfo);
            pt.imData = pt2.imData;
        catch me
            warning('ITP:DTI' , me.getReport);
            return;
        end
    else
        return;
    end
end

x = double(pt.imOffset(1)-1) + (1:size(pt.imData,2));
y = double(pt.imOffset(2)-1) + (1:size(pt.imData,1));

if (~isempty(camcalinfo))
    rpts = camcalinfo.realPtsFromCamPts([x(1) x(end); y(1) y(end)]);
    x = linspace(rpts(1,1), rpts(1,end), length(x));
    y = linspace(rpts(2,1), rpts(2,end), length(y));
    
    
end


if (clipToContour && ~isempty(contour))
    pointsOutside = 2;
    [xx,yy] = meshgrid(x,y);
    mask = reshape(inpolygon(xx(:), yy(:), contour(1,:), contour(2,:)), size(xx));
    mask = imdilate(mask, ones(2*pointsOutside+1));
else
    mask = true(size(pt.imData));
end

cm = gray(64);

if (isempty(underlayImData))
    im = pt.imData;
    im(~mask) = 0;
else
    if (underlayScale > 0)
        underlayImData.im = double(underlayImData.im)/double(max(0.01,max(underlayImData.im(:))))*underlayScale;
        im = pt.overlayImageOnExisting(underlayImData.x, underlayImData.y, zeros(size(underlayImData.im)), camcalinfo, 'mode', 'max', 'clipToContour',clipToContour, 'contour', contour);
        im = im/max(im(:));
        im = max(underlayImData.im, im);
    else
        im = pt.overlayImageOnExisting(underlayImData.x, underlayImData.y, double(underlayImData.im), camcalinfo, 'mode', 'max', 'clipToContour',clipToContour, 'contour', contour);
        %[xx,yy] = meshgrid(x,y);
        %uim = interp2(underlayImData.x, underlayImData.y, underlayImData.im, xx, yy,'*nearest');
        %im = double(pt.imData) + double(uim);
        
        im = im/max(im(:));
    end
    x = underlayImData.x;
    y = underlayImData.y;
end

if (scale == 1 || ~isempty(underlayImData))
    if (size(im,3) > 1)
        set(image (x,y,double(im)), 'Parent', Axes);
    else
        set(imagesc (x,y,double(im)), 'Parent', Axes, patchOptions{:}); shading (Axes, 'flat'); colormap(Axes, cm);
    end
else
    im2 = imresize(double(im), scale);
    x2 = interp1(x, (0:(size(im2,2)- 1))*length(x)/size(im2,2) + 1, 'linear');
    y2 = interp1(y, (0:(size(im2,1)- 1))*length(y)/size(im2,1) + 1, 'linear');
    set(pcolor(Axes, x2,y2,im2), patchOptions{:}); shading (Axes, 'interp'); colormap (Axes, cm);
end

if (showCov)
    hold(Axes, 'on');
    axes(Axes);
    c = [pt.cov(1) pt.cov(2);pt.cov(2) pt.cov(3)];
    [v,d] = eig(c);
    plot (pt.loc(1),pt.loc(2), 'r.', 'MarkerSize', 20)
    set(ellipse(3*sqrt(d(1,1)), 3*sqrt(d(2,2)), atan2(v(2,1),v(1,1)), pt.loc(1), pt.loc(2),'r-'),'LineWidth',3);
    hold(Axes, 'off');
end

if (isempty(cropSize))
    axis (Axes, 'equal')
    axis (Axes, 'tight')
else
    c = [mean(x) mean(y)];
    axis(Axes, [c(1)-cropSize(1)/2, c(1)+cropSize(1)/2,c(2)-cropSize(2)/2, c(2)+cropSize(2)/2]);
end
