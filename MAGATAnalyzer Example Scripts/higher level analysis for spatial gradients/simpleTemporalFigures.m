function axeshandles = simpleTemporalFigures (ad, plot_options, varargin)
%function axeshandles = simpleTemporalFigures (ad, plot_options, varargin)
%
% spatial_navigation_options -- 1 set of options for all esets
% just plots speed vs time, reo rate vs time, and if provided in file,
% temperature vs. time
%

SaveDirectory = [];
forprinting = false;
showlegend = false;
whichGraphs = {};
fignum = 1;
backgroundColor = [];
legendLocation = 'BestOutside';
showtitle = true;
vectorgraphics = [];
bigfont = false;
figureTitle = 'Experiment';

varargin = assignApplicable(varargin);

f = figure(fignum);
clf(f);
ss = get(0, 'ScreenSize');

screenwidth = ss(3); screenheight = ss(4);
figratio = 8/10;
r = 0.8;
figheight = screenheight*r;
figwidth = figheight*figratio;
figpos = round([max(0,(screenwidth - figwidth*1.33)), (screenheight-figheight)/2, figwidth, figheight]);
    
set(f, 'Position', figpos);
get(f, 'Position');  %calling get f,position here prevents a weird bug where text is incorrectly spaced later
% 
% if (~all(get(f,'Position') == figpos))
%     rerun = false;
% else
%     rerun = false;
% end
%get(f, 'Position')
set(f, 'PaperType', 'usletter', 'PaperPosition', [0.25 0.25 8 10], 'PaperOrientation', 'portrait', 'color', 'w', 'inverthardcopy', 'off');
%get(f, 'Position')

leftmargin = 1/8;
rightmargin = .5/8;
allaxeswidth = 1 - leftmargin - rightmargin;
topmargin = 1/11;
bottommargin = 1/11;
h0 = 1-topmargin;
allaxesheight = h0-bottommargin;

wspace2 = 0.2 * allaxeswidth;
wspace3 = 0.12 * allaxeswidth;
wspace4 = 0.12 * allaxeswidth;

hspace = .75/11;

if (isfield(ad, 'temperature_vs_time'))
    nrows = 3;
else
    nrows = 2;
end

h = (allaxesheight - (nrows-1)*hspace)/nrows;
dh = h + hspace;
w1 = allaxeswidth;

w2 = min((allaxeswidth - wspace2)/2, h*1.61803399);
w3 = (allaxeswidth - 2*wspace3)/3;
w4 = (allaxeswidth - 3*wspace4)/4;

centerx = 0.5*(1+leftmargin -rightmargin);

lx2 = centerx - wspace2/2 - w2;
rx2 = centerx + (wspace2)/2;

lx3 = centerx - wspace3 - 3*w3/2;
cx3 = centerx - w3/2;
rx3 = centerx + w3/2 + wspace3;


lx4 = centerx - 3*wspace4/2 - 2*w4;
clx4 = centerx - wspace4/2 - w4;
crx4 = centerx + (wspace4)/2;
rx4 = centerx + 3*(wspace4)/2+w4;

default_properties = get(0,'default');


vectorgraphics = true;

default_properties = get(0,'default');



ccc = 'bgrcymk';
sss = 'sodvh>p^<';
po.lineWidth = 1;
po.color = ccc(1);

po.marker = 'none';
po.plotOptions = {};
po.shadedErrorRegion = true;
po.useGauss = false;
po.directions = [-180,-90, 0, 90];
po.colors = {[0.1 0 1], [0.1 0.8 0.2], [1 0 0], [0.6 0.6 0.0]};
po.labelRadii = false;
po.textOnPolar = true;
po.font = 'Arial';
po.fontsize = 7;

if (nargin == 0)
    axeshandles = po;
    return;
end


spatial_navigation_options = ad.sno;
if (isstruct(spatial_navigation_options))
    fn = fieldnames(spatial_navigation_options);
    for j = 1:length(fn)
        sno.(fn{j}) = spatial_navigation_options.(fn{j});
    end
end

nexp = length(ad);
if (nexp > 1)
    error ('call on only one experiment analyzed data set');
end

po.legendEntry = cellfun(@(x) ['to ' num2str(x) '$^\circ$'], num2cell(po.directions),'UniformOutput',false);
if (nargin > 2 && isstruct(plot_options))
    fn = fieldnames(plot_options);
    for j = 1:length(fn)
        po.(fn{j}) = plot_options.(fn{j});
    end
end
for j = 1:length(po.colors)
    if (ischar(po.colors{j}))
        po.colors{j} = char2rgb(po.colors{j});
    end
end


fontsize = po.fontsize;
font = po.font;
linewidth = po.lineWidth;
set(0,'DefaultAxesFontSize', fontsize);
set(0,'DefaultAxesFontName', font);
set(0, 'DefaultTextInterpreter', 'Tex');
axesopts = {'Box', 'off', 'LineWidth', linewidth, 'FontSize', fontsize, 'FontName', font};

figureTitle(figureTitle == '_') = '-';
annotation ('textbox', [leftmargin, 1-hspace, 1-leftmargin, hspace], 'String', figureTitle, ...
        'FontName', font, 'FontSize', 18, 'Color','k', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top',...
        'LineStyle', 'none');
    
an = 1;
ah(an).name = 'speed_vs_time';
ah(an).pos = [leftmargin, h0-h, w1,h];
ah(an).axes = axes('Position', ah(an).pos, axesopts{:});
shadedErrorPlot(ad.time_axis, ad.speed_vs_time*60, 60*ad.speed_vs_time_eb, [],[0 0 0], 'lineOptions', {'linewidth', 2});
xlabel ('time (s)');
ylabel ('speed (cm/min)');


an = an+1;
ah(an).name = 'reo_vs_time';
ah(an).pos = ah(an-1).pos; ah(an).pos(2) = ah(an).pos(2) - dh;
ah(an).axes = axes('Position', ah(an).pos, axesopts{:});
shadedErrorPlot(ad.time_axis, ad.reo_vs_time, ad.reo_vs_time_eb, [], [ 0 0 0], 'lineOptions', {'linewidth', 2});
xlabel ('time (s)');
ylabel ('reorientation rate (min^{-1})', 'Interpreter', 'Tex');

if (isfield(ad, 'temperature_vs_time'))
    an = an+1;
    ah(an).name = 'temperatureVsTime';
    ah(an).pos = ah(an-1).pos; ah(an).pos(2) = ah(an).pos(2) - dh;
    ah(an).axes = axes('Position', ah(an).pos, axesopts{:});
    plot(ad.time_axis, ad.temperature_vs_time, 'k', 'LineWidth', 2);
    xlabel ('time (s)');
    ylabel ('remperature (^{\circ} C)', 'Interpreter', 'Tex');
end
