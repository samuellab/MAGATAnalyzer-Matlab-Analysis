function [cax, cax2] = polarBackground(radii, cax, locOfZero, varargin)
%function [cax, cax2] = polarBackground(radii, cax, locOfZero, varargin)
%function [cax, cax2] = polarBackground(maxRadius, cax, locOfZero, varargin)
%radii = vector of radii that will be plotted, or the maximum Radius
%if second argument is requested, a new axes is created on top of cax and
%set to be clear with no box
%if locOfZero is passed then zero is rotated to occupy the angle (in
%degrees) specified
%
%varargin:
%'labelRadii' true/[false]:  if true the dotted circles are labeled with
%               their radii
%'notext' [false]/true : if true, no labels for angles or radii
%'numlabels' : number of theta positions to label;
labelRadii = false;
notext = false;
numlabels = 12;
varargin = assignApplicable(varargin);

existsAndDefault('cax', gca);
existsAndDefault('locOfZero', 0);
cax = newplot(cax);

next = lower(get(cax,'NextPlot'));
hold_state = ishold(cax);

% get x-axis text color so grid is in same color
tc = get(cax,'xcolor');
ls = get(cax,'gridlinestyle');

% Hold on to current Text defaults, reset them to the
% Axes' font attributes so tick marks use them.
fAngle  = get(cax, 'DefaultTextFontAngle');
fName   = get(cax, 'DefaultTextFontName');
fSize   = get(cax, 'DefaultTextFontSize');
fWeight = get(cax, 'DefaultTextFontWeight');
fUnits  = get(cax, 'DefaultTextUnits');
set(cax, 'DefaultTextFontAngle',  get(cax, 'FontAngle'), ...
    'DefaultTextFontName',   get(cax, 'FontName'), ...
    'DefaultTextFontSize',   get(cax, 'FontSize'), ...
    'DefaultTextFontWeight', get(cax, 'FontWeight'), ...
    'DefaultTextUnits','data')
rho = radii;
% only do grids if hold is off
if ~hold_state

% make a radial grid
    hold(cax,'on');
% ensure that Inf values don't enter into the limit calculation.
    arho = abs(rho(:));
    maxrho = max(arho(arho ~= Inf));
    mult = 0;
    while (maxrho < 1)
        maxrho = maxrho*10;
        mult = mult+1;
    end
    maxrho = ceil(maxrho*5)/5 / 10^mult;
    hhh=line([-maxrho -maxrho maxrho maxrho],[-maxrho maxrho maxrho -maxrho],'parent',cax);
    axis(cax, 'tight');
    set(cax,'dataaspectratio',[1 1 1],'plotboxaspectratiomode','auto')
    v = [get(cax,'xlim') get(cax,'ylim')];
    ticks = sum(get(cax,'ytick')>=0);
    delete(hhh);
% check radial limits and ticks
    rmin = 0; rmax = v(4); rticks = max(ticks-1,2);
    if rticks > 5   % see if we can reduce the number
        if rem(rticks,2) == 0
            rticks = rticks/2;
        elseif rem(rticks,3) == 0
            rticks = rticks/3;
        end
    end

% define a circle
    th = 0:pi/50:2*pi;
    xunit = cos(th);
    yunit = sin(th);
% now really force points on x/y axes to lie on them exactly
    inds = 1:(length(th)-1)/4:length(th);
    xunit(inds(2:2:4)) = zeros(2,1);
    yunit(inds(1:2:5)) = zeros(3,1);
% plot background if necessary
    if ~ischar(get(cax,'color')),
       patch('xdata',xunit*rmax,'ydata',yunit*rmax, ...
             'edgecolor',tc,'facecolor',get(cax,'color'),...
             'handlevisibility','off','parent',cax);
    end

% draw radial circles
    c82 = cos(82*pi/180);
    s82 = sin(82*pi/180);
    rinc = (rmax-rmin)/rticks;
    for i=(rmin+rinc):rinc:rmax
        hhh = line(xunit*i,yunit*i,'linestyle',ls,'color',tc,'linewidth',1,...
                   'handlevisibility','off','parent',cax);
               
           if (labelRadii && ~notext)
               text((i+rinc/20)*c82,(i+rinc/20)*s82, ...
                   ['  ' num2str(i,2)],'verticalalignment','bottom',...
                   'handlevisibility','off','parent',cax)
           end
    end
    set(hhh,'linestyle','-') % Make outer circle solid

% plot spokes
    th = ((1:6)*2*pi/12) + deg2rad(locOfZero);
    cst = cos(th); snt = sin(th);
    cs = [-cst; cst];
    sn = [-snt; snt];
    line(rmax*cs,rmax*sn,'linestyle',ls,'color',tc,'linewidth',1,...
         'handlevisibility','off','parent',cax)

% annotate spokes in degrees
    rt = 1.05*rmax;
    ddtt = 360/numlabels;
    tthh = (0:(numlabels-1))*ddtt;
    for i = 1:length(tthh)
        if (~notext)
            if (tthh(i) == 180)
                loc = '\pm 180';
            else
                loc = int2str(mod(tthh(i) + 180,360)-180);
            end
            if (cosd(tthh(i)+locOfZero) > 1/sqrt(2))
                ha = 'left';
            else
                if (cosd(tthh(i)+locOfZero) < -1/sqrt(2))
                    ha = 'right';
                else
                    ha = 'center';
                end
            end
            if (sind(tthh(i)+locOfZero) > 1/sqrt(2))
                va = 'bottom';
            else
                if (sind(tthh(i)+locOfZero) < -1/sqrt(2))
                    va = 'top';
                else
                    va = 'middle';
                end
            end
            text(rt*cosd(tthh(i)+locOfZero),rt*sind(tthh(i)+locOfZero),loc,'horizontalalignment',ha,'verticalalignment',va,...
                 'handlevisibility','off','parent',cax)
        end
    end

% set view to 2-D
    view(cax,2);
% set axis limits
    axis(cax,rmax*[-1 1 -1.15 1.15]);
end

% Reset defaults.
set(cax, 'DefaultTextFontAngle', fAngle , ...
    'DefaultTextFontName',   fName , ...
    'DefaultTextFontSize',   fSize, ...
    'DefaultTextFontWeight', fWeight, ...
    'DefaultTextUnits',fUnits );

if ~hold_state
    set(cax,'dataaspectratio',[1 1 1]), axis(cax,'off'); set(cax,'NextPlot',next);
end

if (nargout > 1)
    cax2 = cloneaxes(cax);
    set(cax2, 'Color', 'none', 'XTick', [], 'YTick', [],'dataaspectratio',[1 1 1],'plotboxaspectratiomode','auto'); 
    axis(cax2, 'off');
end
