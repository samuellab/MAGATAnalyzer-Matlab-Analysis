function drawTrackImage (pt, camcalinfo, varargin)
%function drawTrackImage (pt, camcalinfo, varargin)
%@MaggotTrackPoint
%
%passing 'pretty', 'true' opens up the following options
%'scale', s (default 3) -- scales the image to higher resolution using
%imresize
%'contourColor'
%'contourWidth'
%    scale = 3;
%    drawContour = true;
%    drawHeadArrow = true;
%    contourColor = 'r-';
%    spineColor = 'y.';
%    contourWidth = 2;
%    mhColor = 'r';
%    tmColor = 'b-';
%    arrowSize = 2;
%    mhWidth = 2;

if (~exist('camcalinfo', 'var'))
    camcalinfo = [];
end

pretty = false;
Axes = [];
varargin = assignApplicable(varargin);
if (isempty(Axes))
    Axes = gca;
end
%{
if (isempty(camcalinfo))
    h = realPtsToCamera(pt.head, camcalinfo);   
    m = realPtsToCamera(pt.mid, camcalinfo);
    t = realPtsToCamera(pt.tail, camcalinfo);
    c = realPtsToCamera(pt.contour, camcalinfo);
    sp = realPtsToCamera(pt.spine, camcalinfo);
else
  %}
h = pt.head;
m = pt.mid;
t = pt.tail;
c = pt.contour;
sp = pt.spine;
c(:,end+1) = c(:,1); %complete contour

if (~pretty)
    drawTrackImage@ImTrackPoint(pt, camcalinfo, 'Axes', Axes, varargin{:});
    
    hold (Axes, 'on');
  %  plot (Axes, h(1),h(2),'g*', t(1),t(2),'rh', [t(1),m(1),h(1)], [t(2),m(2),h(2)],'y-', c(1,:), c(2,:), 'r-');
    if (~isempty(sp) && all(isfinite(sp(:))))
        plot(Axes, h(1),h(2),'g*', t(1),t(2),'rh', sp(1,:), sp(2,:),'y.-', c(1,:), c(2,:), 'r-');
    else
        plot (Axes, h(1),h(2),'g*', t(1),t(2),'rh', [t(1),m(1),h(1)], [t(2),m(2),h(2)],'y-', c(1,:), c(2,:), 'r-');
    end
    hold (Axes, 'off');
else
    scale = 3;
    drawContour = true;
    drawSpine = true;
    drawHeadArrow = true;
    contourColor = 'r-';
    spineColor = 'y.';
    contourWidth = 2;
    mhColor = 'r';
    tmColor = 'b-';
    spineMarkerSize = 5;
    spineLineWidth = 2;
    arrowSize = 2;
    mhWidth = 2;
    varargin = assignApplicable(varargin);
    drawTrackImage@ImTrackPoint(pt, camcalinfo, 'Axes', Axes, 'scale', scale, varargin{:});
    shading(Axes, 'interp');
    ih = ishold(Axes);
    hold (Axes, 'on');
    if (drawContour)
        c2 = interp1(c', 1:(1/scale):length(c));
        c2 = lowpass1D(c2', scale);        
        plot (Axes, c2(1,[1:end 1]), c2(2,[1:end 1]), contourColor, 'LineWidth', contourWidth);
         
    end
    if (drawSpine && pt.htValid)
        if (~isempty(sp) && all(isfinite(sp(:))))
            plot (Axes, sp(1,:), sp(2,:), spineColor, 'MarkerSize', spineMarkerSize, 'LineWidth', spineLineWidth);
        end
    end
    if (drawHeadArrow && pt.htValid)
        quiver(Axes, m(1), m(2), h(1)-m(1), h(2)-m(2), 0,mhColor,'LineWidth',mhWidth,'MaxHeadSize', arrowSize);
        plot (Axes, [t(1) m(1)], [t(2) m(2)], tmColor, 'LineWidth', mhWidth);
    end
    if (~ih)
        hold (Axes, 'off');
    end
end