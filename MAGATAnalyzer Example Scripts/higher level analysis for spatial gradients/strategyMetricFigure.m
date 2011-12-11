function axeshandles = strategyMetricFigure (ad, plot_options, varargin)
%function axeshandles = strategyMetricFigure (ad, plot_options, varargin)
%
% spatial_navigation_options -- 1 set of options for all esets
% plot_options -- each eset gets its own options
% plot_options.
%    lineWidth -- width of lines
%    color -- color for line and marker
%    legendEntry -- what to put in legend
%    marker -- marker 
%    plotOptions -- additional options to pass to plotting function
%
%   varargin: 'SaveDirectory', directory to save images in (if empty,
%   nothing saved
%   forprinting = false;
%   showlegend = false;
%   whichGraphs = {}; all if empty - choose from {'DirectionHistogram',
%   'RunStartHistogram','ReorientationRateVsHeading',
%   'InstantaneousDeltaThetaVsTheta','SpeedVsDirection',
%   'ReoDirVsHeading','RunDirVsHeading','ReoMagVsHeading',  'HeadSwingAcceptanceHandedness',
%   'FirstHeadSwingHandedness', 'ReoDirDistribution','RunLengthHistogram',
%   'RunDirDistribution','ReoDirDistributionPolar',
%   'RunDirDistributionPolar', 'FirstHeadSwingHandednessPerp'

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

dxc = 0.1;
allaxeswidth = 0.95 - dxc;
allaxesheight = 0.9;
dh = (.975 - allaxesheight);

wspace = 0.1 * allaxeswidth;

w = (allaxeswidth - 2*wspace)/3;
h = 0.18 * allaxesheight;


w4 = w*3/4*(allaxeswidth-dxc)/allaxeswidth;
h4 = w4*8/10;

w4s = (allaxeswidth-dxc-4*w4)/3;


%dw = (1 - (3*w + 2*wspace))/2;

lx = dxc;
cx = dxc + allaxeswidth/2-w/2;
rx = allaxeswidth+dxc-w;

lx4 = 2*dxc;
clx4 = lx4 + w4+w4s;
crx4 = clx4 + w4+w4s;
rx4 =  crx4 + w4+w4s;
x4 = [lx4, clx4, crx4, rx4];

hspace = (allaxesheight - 3*h -2*h4)/4.5;


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
annotation ('textbox', [dxc, 1-dh, 1-dxc, dh], 'String', figureTitle, ...
        'FontName', font, 'FontSize', 18, 'Color','k', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top',...
        'LineStyle', 'none');
    
an = 1;
ah(an).name = 'Rose';
ah(an).pos = [allaxeswidth+dxc-dh, 1-dh*9/10, dh, dh*7/10];
ah(an).axes = axes('Position', ah(an).pos, axesopts{:});
polarBackground(1, ah(an).axes, 0, 'labelRadii', false, 'notext', false, 'numlabels',4);
th = -45:45;
for j = 1:4
    x = [0 cosd(po.directions(j) + th) 0];
    y = [0 sind(po.directions(j) + th) 0];
    patch(x,y,1,'FaceColor', po.colors{j}*0.25 + .75, 'EdgeColor', 'k');
end

an = an+1;
ah(an).name = 'DirectionHistogram';
ah(an).pos = [lx 1-dh-h w h];
ah(an).axes = axes('Position', ah(an).pos, axesopts{:});
quadrantGraph(ad, 'txc', 'thetahist', po, 1/sum(ad.thetahist), true);
ylabel({'Relative probability', 'of orientation'}, 'Interpreter', 'Tex');
xlabel('run heading (^\circ)', 'Interpreter', 'Tex');
set(ah(an).axes, 'XTick', [-180 -90 0 90 180], axesopts{:});


an = an+1;
ah(an).name = 'SpeedVsDirection';
ah(an).pos = ah(an-1).pos; ah(an).pos(1) = cx;
ah(an).axes = axes('Position', ah(an).pos, axesopts{:});

yfield = 'speedVsDir';
ymult = 60;
if (po.useGauss)
    xfield = 'txf';
    yfield = [yfield '_gauss'];
    iscirc = false;
else
    xfield = 'txc';
    iscirc = true;
end
quadrantGraph(ad,xfield,yfield, po, ymult, iscirc);
ylabel ('Run speed (cm/min)', 'Interpreter', 'Tex');
xlabel('run heading (\circ)', 'Interpreter', 'Tex');
set(ah(an).axes, 'XTick', [-180 -90 0 90 180], axesopts{:});

an = an+1;
ah(an).name = 'RunDirVsHeading';
ah(an).pos = ah(an-1).pos; ah(an).pos(1) = rx;
ah(an).axes = axes('Position', ah(an).pos, axesopts{:});
yfield = 'meanrunchange';
ymult = 180/pi;
if (po.useGauss)
    xfield = 'txf';
    yfield = [yfield '_gauss'];
    iscirc = false;
else
    xfield = 'reotxc';
    iscirc = true;
end

quadrantGraph(ad,xfield,yfield, po, ymult, iscirc);
xlabel ('Run start heading (\circ)', 'Interpreter', 'Tex');
ylabel ({'Mean heading change (\circ)', 'in run'}, 'Interpreter', 'Tex');
set(ah(an).axes, 'XTick', [-180 -90 0 90 180], axesopts{:});

an = an+1;
ah(an).name = 'ReorientationRateVsHeading';
ah(an).pos = ah(an-1).pos; ah(an).pos(1) = lx; ah(an).pos(2) = ah(an).pos(2)-h-hspace;
ah(an).axes = axes('Position', ah(an).pos, axesopts{:});

quadrantGraph(ad, 'txc', 'reohist', po, 1, true);
xlabel ('Run heading (\circ)', 'Interpreter', 'Tex');
ylabel ('Turn rate (min^-^1)', 'Interpreter', 'Tex');
set(ah(an).axes, 'XTick', [-180 -90 0 90 180], axesopts{:});


an = an+1;
ah(an).name = 'ReoMagVsHeading';
ah(an).pos = ah(an-1).pos; ah(an).pos(1) = cx;
ah(an).axes = axes('Position', ah(an).pos, axesopts{:});
yfield = 'reoMag';
ymult = 180/pi;
if (po.useGauss)
    xfield = 'txf';
    yfield = [yfield '_gauss'];
    iscirc = false;
else
    xfield = 'reotxc';
    iscirc = true;
end
ad.([yfield '_sqrt']) = sqrt(ad.(yfield));
ad.([yfield '_sqrt_eb']) = 0.5*sqrt(ad.([yfield '_eb']) ./ ad.(yfield));
quadrantGraph(ad,xfield,[yfield '_sqrt'], po, ymult, iscirc);
xlabel ('Prev run heading (\circ)', 'Interpreter', 'Tex');
ylabel ({'Heading change size', '(rms angle, \circ)'}, 'Interpreter', 'Tex');
set(ah(an).axes, 'XTick', [-180 -90 0 90 180], axesopts{:});


an = an+1;
ah(an).name = 'ReoDirVsHeading';
ah(an).pos = ah(an-1).pos; ah(an).pos(1) = rx;
ah(an).axes = axes('Position', ah(an).pos, axesopts{:});
yfield = 'reoDir';
ymult = 180/pi;
if (po.useGauss)
    xfield = 'txf';
    yfield = [yfield '_gauss'];
    iscirc = false;
else
    xfield = 'reotxc';
    iscirc = true;
end

quadrantGraph(ad,xfield,yfield, po, ymult, iscirc);
xlabel ('Prev run heading (\circ)', 'Interpreter', 'Tex');
ylabel ({'Mean heading change (\circ)', 'during reorientation'}, 'Interpreter', 'Tex');
set(ah(an).axes, 'XTick', [-180 -90 0 90 180], axesopts{:});

ind = find(strcmpi('RunDirVsHeading', {ah.name}));
yl1 = get(ah(ind).axes, 'YLim');
yl2 = get(ah(an).axes, 'YLim');
ylm = max(abs([yl1 yl2]));
set([ah([ind an]).axes], 'YLim', [-ylm ylm]);


[~,I] = max(cosd(po.directions - sno.preferredDirection));
tocolor = po.colors{I}/4 + 0.75;
[~,I] = min(cosd(po.directions - sno.preferredDirection));
fromcolor = po.colors{I}/4 + 0.75;

an = an+1;
ah(an).name = 'ReoDirBias';
ah(an).pos = ah(an-1).pos; ah(an).pos(1) = lx; ah(an).pos(2) = ah(an).pos(2)-h-hspace;
ah(an).axes = axes('Position', ah(an).pos, axesopts{:});

bardat = [ad.simpleMetrics.prob_turn_towards 1-ad.simpleMetrics.prob_turn_towards];
ebdat = ad.simpleMetrics.prob_turn_towards_eb*[1 1];
gnames = {['to ' num2str(ad.sno.preferredDirection)], ['to ' num2str(mod(ad.sno.preferredDirection + 180, 360))]};
hhh = barweb_marc(bardat', ebdat', 0.8, gnames);
set(hhh.bars, 'FaceColor', [1 1 1], 'EdgeColor', [0 0 0], 'LineWidth', linewidth);
set(hhh.bars(1), 'FaceColor', tocolor);
set(hhh.bars(2), 'FaceColor', fromcolor);
set(hhh.ax, 'FontSize', fontsize, 'FontName', font, 'YLim', [0 1], 'Box', 'off', 'LineWidth', linewidth);
ylabel ({'Distribution of turns', 'from perp. direction'});



an = an+1;
ah(an).name = 'FirstHSDirBias';
ah(an).pos = ah(an-1).pos; ah(an).pos(1) = cx; 
ah(an).axes = axes('Position', ah(an).pos, axesopts{:});

bardat = [ad.firstHSProbTowards 1-ad.firstHSProbTowards];
ebdat = ad.firstHSProbTowards_eb*[1 1];
gnames = {['to ' num2str(ad.sno.preferredDirection)], ['to ' num2str(mod(ad.sno.preferredDirection + 180, 360))]};
hhh = barweb_marc(bardat', ebdat', 0.8, gnames);
set(hhh.bars, 'FaceColor', [1 1 1], 'EdgeColor', [0 0 0], 'LineWidth', linewidth);
set(hhh.bars(1), 'FaceColor', tocolor);
set(hhh.bars(2), 'FaceColor', fromcolor);
set(hhh.ax, 'FontSize', fontsize, 'FontName', font, 'YLim', [0 1], 'Box', 'off', 'LineWidth', linewidth);
ylabel ({'Distribution of first head sweeps', 'from perp. direction'});


an = an+1;
ah(an).name = 'HSAcceptanceBias';
ah(an).pos = ah(an-1).pos; ah(an).pos(1) = rx; 
ah(an).axes = axes('Position', ah(an).pos, axesopts{:});

bardat = [ad.headSwingAcceptanceRateTowards ad.headSwingAcceptanceRateAway];
ebdat = [ad.headSwingAcceptanceRateTowards_eb ad.headSwingAcceptanceRateAway_eb];
gnames = {['to ' num2str(ad.sno.preferredDirection)], ['to ' num2str(mod(ad.sno.preferredDirection + 180, 360))]};
hhh = barweb_marc(bardat', ebdat', 0.8, gnames);
set(hhh.bars, 'FaceColor', [1 1 1], 'EdgeColor', [0 0 0], 'LineWidth', linewidth);
set(hhh.bars(1), 'FaceColor', tocolor);
set(hhh.bars(2), 'FaceColor', fromcolor);
set(hhh.ax, 'FontSize', fontsize, 'FontName', font, 'YLim', [0 1], 'Box', 'off', 'LineWidth', linewidth);
ylabel ({'Probability of starting run', 'during head sweep'});



an = an+1;
ah(an).name = 'RunDirDistribution';

pos = [lx4 ah(an-1).pos(2)-h4-hspace*2 w4 h4];
for j = 1:4
    denom = sum(ad.run_dtheta_dist{j});
    maxr(j) = max(ad.run_dtheta_dist{j}+ad.run_dtheta_dist_eb{j})/denom;
end
maxr = max(maxr);
xp = [0 0 .5 .5];
yp = [0.5 0 0.5 0];

for j = 1:4
    pos(1) = x4(j);
    ah(an).axes(j) = axes('Position', pos, axesopts{:});
    [~,ind(j)] = min(abs(mod(po.directions - ad.reobasedirections(j) + 180, 360)- 180));
    denom = sum(ad.run_dtheta_dist{j});% + ad.reo_dtheta_dist_nohs{j});
    bc = [1 1 1];
    fc2 = 0.4*po.colors{ind(j)} + 0.6*bc;
    
    [~,ax2(j)] = polarBackground(maxr, ah(an).axes(j), ad.reobasedirections(j), 'labelRadii', po.labelRadii, 'notext', ~po.textOnPolar, 'numlabels', 4);
    
    
    [hh,heb] = polarBarPlotWError(ad.dtx, ad.run_dtheta_dist{j}/denom, ad.run_dtheta_dist_eb{j}/denom,[], fc2, 'locOfZero', ad.reobasedirections(j), 'FaceColor', fc2, 'EdgeColor', po.colors{ind(j)}, 'Parent', ax2(j), 'LineWidth', 1, 'curveEB', true); 
    set(heb, 'LineWidth', po.lineWidth);
    
    [x,y] = dsxy2figxy_marc(ah(an).axes(j), maxr*cosd(ad.reobasedirections(j))*[-1 -.5],  maxr*sind(ad.reobasedirections(j))*[-1 -.5]);
    annotation('arrow', x,y, 'Color', po.colors{ind(j)}, 'LineWidth', po.lineWidth*1.5);
    
    pp = get(ah(an).axes(j), 'OuterPosition');
    if (cosd(ad.reobasedirections) == -1)
        str = '\pm180\circ';
    else
        str = [num2str(ad.reobasedirections(j)) '\circ'];
    end
    annotation ('textbox', [pos(1)+w4/2, pos(2)+h4+hspace/2, 0.001, 0.001], 'String', str, ...
        'FontName', font, 'FontSize', fontsize+2, 'Color', po.colors{ind(j)}, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom',...
        'LineStyle', 'none');
end

annotation ('textbox', [0, pos(2)+h4+hspace/2, x4(1)-dxc/2, 0.001], 'String', 'Initial heading:', ...
        'FontName', font, 'FontSize', fontsize, 'Color','k', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom',...
        'LineStyle', 'none');

annotation ('textbox', [0, pos(2)+h4/2, x4(1)-dxc/2, 0.001], 'String', {'Heading changes','during runs'}, ...
        'FontName', font, 'FontSize', fontsize, 'Color','k', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'middle',...
        'LineStyle', 'none');
    
    
ah(an).pos = pos;




an = an+1;
ah(an).name = 'ReoDirDistribution';

pos = [lx4 ah(an-1).pos(2)-h4-hspace/2 w4 h4];
for j = 1:4
    denom = sum(ad.reo_dtheta_dist_hs{j});
    maxr(j) = max(ad.reo_dtheta_dist_hs{j}+ad.reo_dtheta_dist_hs_eb{j})/denom;
end
maxr = max(maxr);
xp = [0 0 .5 .5];
yp = [0.5 0 0.5 0];

for j = 1:4
    pos(1) = x4(j);
    ah(an).axes(j) = axes('Position', pos, axesopts{:});
    [~,ind(j)] = min(abs(mod(po.directions - ad.reobasedirections(j) + 180, 360)- 180));
    denom = sum(ad.reo_dtheta_dist_hs{j});% + ad.reo_dtheta_dist_nohs{j});
    bc = [1 1 1];
    fc2 = 0.4*po.colors{ind(j)} + 0.6*bc;
    
    [~,ax2(j)] = polarBackground(maxr, ah(an).axes(j), ad.reobasedirections(j), 'labelRadii', po.labelRadii, 'notext', ~po.textOnPolar, 'numlabels', 4);
    
    
    [hh,heb] = polarBarPlotWError(ad.dtx, ad.reo_dtheta_dist_hs{j}/denom, ad.reo_dtheta_dist_hs_eb{j}/denom,[], fc2, 'locOfZero', ad.reobasedirections(j), 'FaceColor', fc2, 'EdgeColor', po.colors{ind(j)}, 'Parent', ax2(j), 'LineWidth', 1, 'curveEB', true); 
    set(heb, 'LineWidth', po.lineWidth);
    
    [x,y] = dsxy2figxy_marc(ah(an).axes(j), maxr*cosd(ad.reobasedirections(j))*[-1 -.5],  maxr*sind(ad.reobasedirections(j))*[-1 -.5]);
    annotation('arrow', x,y, 'Color', po.colors{ind(j)}, 'LineWidth', po.lineWidth*1.5);
    
end
annotation ('textbox', [0, pos(2)+h4/2, x4(1)-dxc/2, 0.001], 'String', {'Heading changes','during turns'}, ...
        'FontName', font, 'FontSize', fontsize, 'Color','k', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'middle',...
        'LineStyle', 'none');


ah(an).pos = pos;
if (nargout > 0)
    axeshandles = ah;
end
set(0, default_properties);






    
function quadrantGraph(ad, xfield, yfield,po, mult, iscirc)
%function quadrantGraph(ad, xfield, yfield,po, mult, iscirc)
[directions,I] = sort(po.directions);
colors = po.colors(I);
existsAndDefault('mult', 1);

tx = ad.(xfield);


q = zeros(size(tx));
for j = 1:length(tx)
    [~,I] = max(cosd(directions - tx(j)));
    q(j) = I;
end
yd = mult*ad.(yfield);
ebd = mult*ad.([yfield '_eb']);

if (iscirc)
    yd = yd([1:end 1]);
    ebd = ebd([1:end 1]);
end

for j = 1:length(directions)
    inds = find(q == j);
    [xdata{j},I] = sort(tx(inds));
    inds = inds(I);
    ydata{j} = yd(inds);
    ebdata{j} = ebd(inds);
end
for j = 1:length(xdata)
    if (any(diff(xdata{j}) > min(diff(directions))))        
        [~,I] = max(diff(xdata{j}));
        ni = length(xdata) + 1;
        colors{ni} = colors{j};
        xdata{ni} = xdata{j}((I+1):end);
        xdata{j} = xdata{j}(1:I);        
        ydata{ni} = ydata{j}((I+1):end);
        ydata{j} = ydata{j}(1:I);
        ebdata{ni} = ebdata{j}((I+1):end);
        ebdata{j} = ebdata{j}(1:I);
    end
end



xf = zeros(size(xdata)); xr = xf; yf = xf; yr = xf; ebf = xf; ebr = xf;

for j = 1:length(xdata)
    ahead = mod(j, length(xdata)) + 1;
    behind = mod(j-2, length(xdata)) + 1;
   
    xf(j) = 0.5* (xdata{j}(end) + mod(xdata{ahead}(1) - xdata{j}(end), 360) + xdata{j}(end));
    xr(j) = 0.5* (-mod(xdata{j}(1) - xdata{behind}(end), 360) + xdata{j}(1) + xdata{j}(1));
    yf(j) = 0.5* (ydata{j}(end) + ydata{ahead}(1));
    yr(j) = 0.5* (ydata{behind}(end) + ydata{j}(1));
    ebf(j) = 0.5* (ebdata{j}(end) + ebdata{ahead}(1));
    ebr(j) = 0.5* (ebdata{behind}(end) + ebdata{j}(1));
end

if (po.shadedErrorRegion)
    for j = 1:length(xdata)
        xxdata{j} = [xr(j) xdata{j} xf(j)];
        yydata{j} = [yr(j) ydata{j} yf(j)];
        eebdata{j} = [ebr(j) ebdata{j} ebf(j)];
    end
    hh = shadedErrorPlot(xxdata, yydata, eebdata, [], colors, 'LineWidth', po.lineWidth);
    set(hh((length(xdata)+1):end), 'LineWidth', po.lineWidth, po.plotOptions{:});
    hold on;
    if (length(tx) > 60)
        marker = 'none';
    else
        marker = po.marker;
    end
    for j = 1:length(xdata)
        plot (xdata{j}, ydata{j}, 'Color', colors{j},'Marker', marker, 'LineStyle', 'none', 'LineWidth', po.lineWidth, po.plotOptions{:}); 
    end
    hold off;
else
    for j = 1:length(xdata)
        plot ([xr(j) xdata{j} xf(j)], [yr(j) ydata{j} yf(j)], 'k-', 'Color', colors{j}, 'LineWidth', po.lineWidth, po.plotOptions{:}); hold on;
        errorbar (xdata{j}, ydata{j}, ebdata{j}, 'Color', colors{j}, 'Marker', po.marker, 'LineWidth', po.lineWidth, po.plotOptions{:});
    end
    hold off
end

%     
% function quadrantGraph(ad, fieldname,po, mult)
% %function quadrantGraph(ad, fieldname,po, mult)
% directions = po.directions;
% colors = po.colors;
% existsAndDefault('mult', 1);
% q = zeros(size(ad.txc));
% for j = 1:length(ad.txc)
%     [~,I] = max(cosd(directions).*cosd(ad.txc(j)) + sind(directions).*sind(ad.txc(j)));
%     q(j) = I;
% end
% 
% yd = mult*ad.(fieldname)([1:end 1]);
% ebd = mult*ad.([fieldname '_eb'])([1:end 1]);
% 
% 
% for j = 1:(length(ad.txc) - 1)
%     xdata = interp1(ad.txc, [j j+0.5]);
%     ydata = interp1(yd, [j j+0.5]);
%     plot (xdata, ydata, 'k-', 'Color',colors{q(j)}, 'LineWidth', po.lineWidth,  po.plotOptions{:}); hold on;
%     xdata = interp1(ad.txc, [j+0.5 j+1]);
%     ydata = interp1(yd, [j+0.5 j+1]);
%     plot (xdata, ydata, 'k-', 'Color',colors{q(j+1)}, 'LineWidth', po.lineWidth,  po.plotOptions{:}); hold on;
% end
% 
% for j = 1:length(directions)
%     inds = q == j;
%     
%     hhhh(j) = errorbar (ad.txc(inds), yd(inds), ebd(inds), 'k.', 'Color',colors{j}); hold on;
%    % get(hhhh)
%    % bob = get(get(hhhh,'Children'))
% %    pause
% end
% set(hhhh, 'Marker', po.marker, 'LineWidth', po.lineWidth,  po.plotOptions{:});
% 
% function quadrantGraphReo(ad, fieldname,po, mult)
% %function quadrantGraph(ad, fieldname,po, mult)
% directions = po.directions;
% colors = po.colors;
% existsAndDefault('mult', 1);
% q = zeros(size(ad.reotxc));
% for j = 1:length(ad.reotxc)
%     [~,I] = max(cosd(directions).*cosd(ad.reotxc(j)) + sind(directions).*sind(ad.reotxc(j)));
%     q(j) = I;
% end
% 
% yd = mult*ad.(fieldname)([1:end 1]);
% ebd = mult*ad.([fieldname '_eb'])([1:end 1]);
% 
% for j = 1:length(directions)
%     inds = q == j;
%     errorbar (ad.reotxc(inds), yd(inds), ebd(inds), 'k.', 'Color',colors{j}, 'Marker', po.marker, 'LineWidth', po.lineWidth,  po.plotOptions{:}); hold on;
% end
% 
% for j = 1:(length(ad.reotxc) - 1)
%     xdata = interp1(ad.reotxc, [j j+0.5]);
%     ydata = interp1(yd, [j j+0.5]);
%     plot (xdata, ydata, 'k-', 'Color',colors{q(j)}, 'LineWidth', po.lineWidth,  po.plotOptions{:}); hold on;
%     xdata = interp1(ad.reotxc, [j+0.5 j+1]);
%     ydata = interp1(yd, [j+0.5 j+1]);
%     plot (xdata, ydata, 'k-', 'Color',colors{q(j+1)}, 'LineWidth', po.lineWidth,  po.plotOptions{:}); hold on;
% end