function [cax, cax2] = polarBackgroundHalfMoon(radii, cax, locOfZero, centerOfFace1, colorOfFace1, colorOfFace2, varargin)
%function [cax, cax2] = polarBackground(radii, cax, locOfZero)
%function [cax, cax2] = polarBackground(maxRadius, cax, locOfZero)
%radii = vector of radii that will be plotted, or the maximum Radius
%if second argument is requested, a new axes is created on top of cax and
%set to be clear with no box
%if locOfZero is passed then zero is rotated to occupy the angle (in
%degrees) specified
labelRadii = false;
notext = false;
numlabels = 12;
linewidth = 1;
varargin = assignApplicable(varargin);

existsAndDefault('cax', gca);
existsAndDefault('locOfZero', 0);
existsAndDefault('centerOfFace1', 0);
existsAndDefault('colorOfFace1', [1 1 1]);
if (ischar(colorOfFace1))
    colorOfFace1 = char2rgb(colorOfFace1);
end
existsAndDefault('colorOfFace2', 1 - colorOfFace1);
if (ischar(colorOfFace2))
    colorOfFace2 = char2rgb(colorOfFace2);
end
cax = newplot(cax);

next = lower(get(cax,'NextPlot'));
hold_state = ishold(cax);

% get x-axis text color so grid is in same color
tc = get(cax,'xcolor');
%ls = get(cax,'gridlinestyle');
ls = ':';

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
    fc = colorOfFace1;
    tc =  1 - colorOfFace1;
    th = linspace(-pi/2, pi/2, 50) + deg2rad(centerOfFace1);
    xunit = cos(th);
    yunit = sin(th);
    
    patch('xdata',xunit([1:end 1])*rmax,'ydata',yunit([1:end 1])*rmax, ...
        'edgecolor','none','facecolor',fc,...
        'handlevisibility','off','parent',cax);
    
    % draw radial circles
    rang = 82 - mod(locOfZero, 30);
    c82 = cos(rang*pi/180);
    s82 = sin(rang*pi/180);
    rinc = (rmax-rmin)/rticks;
    for i=(rmin+rinc):rinc:rmax
        hhh = line(xunit*i,yunit*i,'linestyle',ls,'color',tc,'linewidth',linewidth,...
            'handlevisibility','off','parent',cax);
        if (labelRadii && ~notext)
            if (cosd(centerOfFace1 - rang) > 0)
                text((i+rinc/20)*c82,(i+rinc/20)*s82, ...'
                    ['  ' num2str(i)],'verticalalignment','bottom',...
                    'handlevisibility','off','parent',cax, 'color', tc)
            end
        end
    end
    % Make outer circle solid
    set(hhh,'linestyle','-', 'color', [0 0 0]);
    
    
    fc = colorOfFace2;
    tc = 1 - colorOfFace2;
    th = pi + linspace(-pi/2, pi/2, 50) + deg2rad(centerOfFace1);
    xunit = cos(th);
    yunit = sin(th);
    
    patch('xdata',xunit([1:end 1])*rmax,'ydata',yunit([1:end 1])*rmax, ...
        'edgecolor','none','facecolor',fc,...
        'handlevisibility','off','parent',cax);
    
    % draw radial circles
    for i=(rmin+rinc):rinc:rmax
        hhh = line(xunit*i,yunit*i,'linestyle',ls,'color',tc,'linewidth',linewidth,...
            'handlevisibility','off','parent',cax);
         if (labelRadii && ~notext)
              
            if (cosd(centerOfFace1 - rang) < 0)
                text((i+rinc/20)*c82,(i+rinc/20)*s82, ...'
                    ['  ' num2str(i)],'verticalalignment','bottom',...
                    'handlevisibility','off','parent',cax, 'color', tc)
            end
         end
    end
    % Make outer circle solid
    set(hhh,'linestyle','-', 'color', [0 0 0]);
    
    
    % plot spokes
    
    for th = ((1:12)*2*pi/12) + deg2rad(locOfZero);
        cst = cos(th); snt = sin(th);
        cs = [0; cst];
        sn = [0; snt];
        if (cos(th - deg2rad(centerOfFace1)) > 0)
            tc = 1 - colorOfFace1;
        else
            tc = 1 - colorOfFace2;
        end
        line(rmax*cs,rmax*sn,'linestyle',ls,'color',tc,'linewidth',linewidth,...
            'handlevisibility','off','parent',cax)
    end
%     %annotate spokes in degrees
%     rt = 1.1*rmax;
%     th = ((1:6)*2*pi/12) + deg2rad(locOfZero);
%     cst = cos(th); snt = sin(th);
%     for i = 1:length(th)
%         if (~notext)
%               
%             text(rt*cst(i),rt*snt(i),int2str(mod(i*30 + 180, 360)-180),...
%                 'horizontalalignment','center',...
%                 'handlevisibility','off','parent',cax);
%         
%             loc = int2str(mod(i*30, 360) - 180);
%             text(-rt*cst(i),-rt*snt(i),loc,'horizontalalignment','center',...
%                 'handlevisibility','off','parent',cax)
%         end
%     end
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
