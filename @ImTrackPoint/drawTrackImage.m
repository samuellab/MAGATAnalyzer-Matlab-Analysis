function drawTrackImage (pt, camcalinfo, varargin)
%function drawTrackImage (pt, camcalinfo, varargin)
%@ImTrackPoint
%
%optional arguments, with defaults
%fid = [];
%scale = 1;
%cropSize = [];
%Axes = [];
%underlayImData = []; underlayImData.x, underlayImData.y, underlayImData.im
%showCov = false;
%
fid = [];
existsAndDefault('camcalinfo', []);
scale = 1;
cropSize = [];
Axes = [];
underlayImData = [];
showCov = false;
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
    rpts = camcalinfo.realPtsFromCamPts([x;y]);
    x = rpts(1,:);
    y = rpts(2,:);
end
    
cm = gray(64);

if (isempty(underlayImData))
    im = pt.imData;
else
    [xx,yy] = meshgrid(x,y);
    uim = interp2(underlayImData.x, underlayImData.y, underlayImData.im, xx, yy,'*nearest');
    im = double(pt.imData) + double(uim);
end

if (scale == 1)
    pcolor (Axes, x,y,double(im)); shading (Axes, 'flat'); colormap(Axes, cm);
else
    im2 = imresize(double(im), scale);
    x2 = interp1(x, (0:(size(im2,2)- 1))*length(x)/size(im2,2) + 1, 'linear');
    y2 = interp1(y, (0:(size(im2,1)- 1))*length(y)/size(im2,1) + 1, 'linear');
    pcolor(Axes, x2,y2,im2); shading (Axes, 'interp'); colormap (Axes, cm);
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
