function simpleMetricFigure(ad, expname, baselinead, baselinename, varargin)
%function simpleMetricFigure(ad, expname, baselinead, baselinename, varargin)


fignum = 1;

font = 'Arial';
fontsize = 9;
bigfontsize = 16;
navindylim = 0.5;
runfracylim = 0.5;
runlenylim = 1.5;
runtimeylim = 40;
speedylim = 5;
lengthylim = 4;
numhsylim = 3;
linewidth = 1;

varargin = assignApplicable(varargin);


hasbaseline = existsAndDefault('baselinead', []);
existsAndDefault('baselinename', 'baseline');
f = figure(fignum);
clf(f);
ss = get(0, 'ScreenSize');

screenwidth = ss(3); screenheight = ss(4);
figratio = 8/10;
r = 0.8;
figheight = screenheight*r;
figwidth = figheight*figratio;
figpos = round([max((screenwidth - 2*figwidth)/3,0), (screenheight-figheight)/2, figwidth, figheight]);
    
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

allaxeswidth = 0.92;
allaxesheight = 0.99;
dw = (1 - allaxeswidth)/2;
dh = (1 - allaxesheight)/2;
spacing = 0.04*allaxesheight;

w = 1 * allaxeswidth;
h = 0.1 * allaxesheight;

textaxespos = [dw 1-dh-h w h];
textaxes = axes('Position', textaxespos, 'Box', 'Off', 'XLim', [0 1], 'YLim', [-1 0], 'XTick', [], 'YTick', []);
axis(textaxes, 'off');

existsAndDefault('expname', 'Experiment Summary');
expname(expname == '_') = '-';
if (hasbaseline)
    expname = [expname ' {\color{blue}(baseline = ' baselinename ')}'];
end
t1 = text(0,0, expname, 'FontName', font, 'FontSize', bigfontsize+2, 'FontWeight', 'bold', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top', 'Parent', textaxes, 'Interpreter', 'Tex');
p = get(t1, 'Extent');
t2 = text(0,p(2), ad.summary_message, 'FontName', font, 'FontSize', bigfontsize-2, 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top', 'Parent', textaxes);

w = 0.45 * allaxeswidth;
h = 0.2 * allaxesheight;
bardat = [ad.navind ad.navind_expt];
ebdat = [ad.navind_eb ad.navind_expt_eb];

gnames = ['all' cellfun(@(x) ['exp ' num2str(x)], num2cell(1:length(ad.navind_expt)), 'UniformOutput', false)];
if (hasbaseline)
    bardat = [bardat baselinead.navind];
    ebdat = [ebdat baselinead.navind_eb];
    gnames = [gnames 'baseline'];
end

parnavaxespos = [dw textaxespos(2)-h-spacing w h];
parnavaxes = axes('Position', parnavaxespos, 'Box', 'Off','FontName', font, 'FontSize', fontsize, 'LineWidth', linewidth);
hhh = barweb_marc(bardat(1,:)', ebdat(1,:)', 0.8, gnames);
%th = rotateticklabel(hhh.ax, 20);
%get(hhh.ax)
set(hhh.ax, 'FontSize', fontsize, 'FontName', font, 'YLim', [-navindylim navindylim], 'Box', 'off', 'LineWidth', linewidth);
if (length(bardat) > 6)
    if (hasebaseline)
        set(hhh.ax, 'XTick', [1 length(bardat)/2 length(bardat)], 'XTickLabel', {'all', 'individual experiments', 'baseline'});
    else
        set(hhh.ax, 'XTick', [1 length(bardat)/2], 'XTickLabel', {'all', 'individual experiments'});
    end
end

set(hhh.bars(1), 'FaceColor', [0.4 0.4 0.4], 'EdgeColor', [0 0 0], 'LineWidth', linewidth);
set(hhh.bars(2:end), 'FaceColor', [0.7 0.7 0.7], 'EdgeColor', [0 0 0], 'LineWidth', linewidth);
if (hasbaseline)
    set(hhh.bars(end), 'FaceColor', [0 0 1], 'EdgeColor', [0 0 0], 'LineWidth', linewidth);
end

xl = get(parnavaxes, 'Xlim'); yl = get(parnavaxes, 'YLim');
text (mean(xl), max(yl), ['Navigation Index towards ' num2str(ad.sno.preferredDirection) '$^{\circ}$'], 'Interpreter', 'Latex', 'FontName', font, 'FontSize', fontsize, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'FontWeight', 'bold');



perpnavaxespos = [1-dw-w textaxespos(2)-h-spacing w h];
perpnavaxes = axes('Position', perpnavaxespos, 'Box', 'Off','FontName', font, 'FontSize', fontsize, 'LineWidth', linewidth);
hhh = barweb_marc(bardat(2,:)', ebdat(2,:)', 0.8, gnames);
%th = rotateticklabel(hhh.ax, 20);
%get(hhh.ax)
set(hhh.ax, 'FontSize', fontsize, 'FontName', font, 'YLim', [-navindylim navindylim], 'Box', 'off', 'LineWidth', linewidth);
if (length(bardat) > 6)
    if (hasebaseline)
        set(hhh.ax, 'XTick', [1 length(bardat)/2 length(bardat)], 'XTickLabel', {'all', 'individual experiments', 'baseline'});
    else
        set(hhh.ax, 'XTick', [1 length(bardat)/2], 'XTickLabel', {'all', 'individual experiments'});
    end
end
set(hhh.bars(1), 'FaceColor', [0.4 0.4 0.4], 'EdgeColor', [0 0 0], 'LineWidth', linewidth);
set(hhh.bars(2:end), 'FaceColor', [0.7 0.7 0.7], 'EdgeColor', [0 0 0], 'LineWidth', linewidth);
if (hasbaseline)
    set(hhh.bars(end), 'FaceColor', [0 0 1], 'EdgeColor', [0 0 0], 'LineWidth', linewidth);
end

xl = get(perpnavaxes, 'Xlim'); yl = get(perpnavaxes, 'YLim');
text (mean(xl), max(yl), ['Navigation Index towards ' num2str(mod(ad.sno.preferredDirection + 90, 360)) '$^{\circ}$'], 'Interpreter', 'Latex', 'FontName', font, 'FontSize', fontsize, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'FontWeight', 'bold');



w = 0.27*allaxeswidth;
h = 0.18*allaxeswidth;

sm = ad.simpleMetrics;
gnames = cellfun(@(x) num2str(x), num2cell(sm.quadrants), 'UniformOutput', false);


pos = [dw perpnavaxespos(2)-h-spacing w h];
nrunsaxes = axes('Position', pos, 'Box', 'Off', 'FontName', font, 'FontSize', fontsize, 'LineWidth', linewidth);
bardat = sm.numRuns/sum(sm.numRuns);
ebdat = sqrt(bardat.*(1-bardat)./sm.numRuns);
hhh = barweb_marc(bardat', ebdat', 0.8, gnames);
set(hhh.bars, 'FaceColor', [1 1 1], 'EdgeColor', [0 0 0], 'LineWidth', linewidth);
set(hhh.ax, 'FontSize', fontsize, 'FontName', font, 'YLim', [0 runfracylim], 'Box', 'off', 'LineWidth', linewidth);
xl = get(nrunsaxes, 'Xlim'); yl = get(nrunsaxes, 'YLim');
if (hasbaseline)
    hold(nrunsaxes, 'on')
    bardat = baselinead.simpleMetrics.numRuns/sum(baselinead.simpleMetrics.numRuns);
    ebdat = sqrt(bardat.*(1-bardat)./baselinead.simpleMetrics.numRuns);
    errorbar(nrunsaxes, 1:4, bardat, ebdat, 'b*', 'LineWidth', linewidth);
    hold(nrunsaxes, 'off');
end
text (mean(xl), max(yl), 'fraction of runs / quadrant', 'FontName', font, 'FontSize', fontsize, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');


pos = [0.5-w/2 perpnavaxespos(2)-h-spacing w h];
runlenaxes = axes('Position', pos, 'Box', 'Off', 'FontName', font, 'FontSize', fontsize, 'LineWidth', linewidth);
bardat = sm.runLen;
ebdat = sm.runLen_eb;
hhh = barweb_marc(bardat', ebdat', 0.8, gnames);
set(hhh.bars, 'FaceColor', [1 1 1], 'EdgeColor', [0 0 0], 'LineWidth', linewidth);
set(hhh.ax, 'FontSize', fontsize, 'FontName', font, 'YLim', [0 runlenylim], 'Box', 'off', 'LineWidth', linewidth);
xl = get(runlenaxes, 'Xlim'); yl = get(runlenaxes, 'YLim');
t1 = text (mean(xl), max(yl), 'mean run length by quadrant', 'FontName', font, 'FontSize', fontsize, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
p = get(t1, 'Extent');
text (mean(xl), p(2), '(cm)', 'FontName', font, 'FontSize', fontsize, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top');
if (hasbaseline)
    hold(runlenaxes, 'on')
    bardat = baselinead.simpleMetrics.runLen;
    ebdat = baselinead.simpleMetrics.runLen_eb;
    errorbar(runlenaxes, 1:4, bardat, ebdat, 'b*', 'LineWidth', linewidth);
    hold (runlenaxes, 'off');
end

pos = [1-dw-w perpnavaxespos(2)-h-spacing w h];
rtpos = pos;
runtimaxes = axes('Position', pos, 'Box', 'Off', 'FontName', font, 'FontSize', fontsize, 'LineWidth', linewidth);
bardat = sm.runTime;
ebdat = sm.runTime_eb;
hhh = barweb_marc(bardat', ebdat', 0.8, gnames);
set(hhh.bars, 'FaceColor', [1 1 1], 'EdgeColor', [0 0 0], 'LineWidth', linewidth);
set(hhh.ax, 'FontSize', fontsize, 'FontName', font, 'YLim', [0 runtimeylim], 'Box', 'off', 'LineWidth', linewidth);
xl = get(runtimaxes, 'Xlim'); yl = get(runtimaxes, 'YLim');
t1 = text (mean(xl), max(yl), 'mean run time by quadrant', 'FontName', font, 'FontSize', fontsize, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
p = get(t1, 'Extent');
text (mean(xl), p(2), '(s)', 'FontName', font, 'FontSize', fontsize, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top');
if (hasbaseline)
    hold(runtimaxes, 'on')
    bardat = baselinead.simpleMetrics.runTime;
    ebdat = baselinead.simpleMetrics.runTime_eb;
    errorbar(runtimaxes, 1:4, bardat, ebdat, 'b*', 'LineWidth', linewidth);
    hold (runtimaxes, 'off');
end


pos = [dw rtpos(2)-h-spacing w h];
reobiasaxes = axes('Position', pos, 'Box', 'Off', 'FontName', font, 'FontSize', fontsize, 'LineWidth', linewidth);
bardat = [sm.prob_turn_towards 1-sm.prob_turn_towards];
ebdat = sm.prob_turn_towards_eb*[1 1];
gnames = {['to ' num2str(ad.sno.preferredDirection)], ['to ' num2str(mod(ad.sno.preferredDirection + 180, 360))]};
hhh = barweb_marc(bardat', ebdat', 0.8, gnames);
set(hhh.bars, 'FaceColor', [1 1 1], 'EdgeColor', [0 0 0], 'LineWidth', linewidth);
set(hhh.ax, 'FontSize', fontsize, 'FontName', font, 'YLim', [0 .8], 'Box', 'off', 'LineWidth', linewidth);
if (hasbaseline)
    hold(reobiasaxes, 'on')
    bardat = [baselinead.simpleMetrics.prob_turn_towards 1-baselinead.simpleMetrics.prob_turn_towards];
    ebdat = baselinead.simpleMetrics.prob_turn_towards_eb*[1 1];
    errorbar(reobiasaxes, bardat, ebdat, 'b*', 'LineWidth', linewidth);
    hold (reobiasaxes, 'off');
end
xl = get(reobiasaxes, 'Xlim'); yl = get(reobiasaxes, 'YLim');
t1 = text (mean(xl), max(yl), 'prob turn from perp', 'FontName', font, 'FontSize', fontsize, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
p = get(t1, 'Extent');
text (mean(xl), p(2), 'is to direction', 'FontName', font, 'FontSize', fontsize, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top');

pos = [0.5-w/2 rtpos(2)-h-spacing w h];
runspeedaxes = axes('Position', pos, 'Box', 'Off', 'FontName', font, 'FontSize', fontsize, 'LineWidth', linewidth);
bardat = [sm.meanspeed sm.meanrunspeed]*60;
hhh = bar(1:2, bardat, 0.8, 'FaceColor', [1 1 1], 'EdgeColor', [0 0 0], 'LineWidth', linewidth);
set(runspeedaxes, 'FontSize', fontsize, 'FontName', font, 'YLim', [0 speedylim], 'Box', 'off','XTick', [1 2], 'XTickLabel', {'all time', 'in runs'}, 'LineWidth', linewidth);
xl = get(runspeedaxes, 'Xlim'); yl = get(runspeedaxes, 'YLim');
t1 = text (mean(xl), max(yl), 'mean speed', 'FontName', font, 'FontSize', fontsize, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
p = get(t1, 'Extent');
text (mean(xl), p(2), '(cm/min)', 'FontName', font, 'FontSize', fontsize, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top');
if (hasbaseline)
    hold(runspeedaxes, 'on')
    bardat = [baselinead.simpleMetrics.meanspeed baselinead.simpleMetrics.meanrunspeed]*60;
    plot(runspeedaxes, 1:2, bardat,'b*', 'LineWidth', linewidth);
    hold (runspeedaxes, 'off');
end

if (hasbaseline)
    pos = [1-dw-w rtpos(2)-h-spacing w h];
    rmpos = pos;
    fracrunsaxes = axes('Position', pos, 'Box', 'Off', 'FontName', font, 'FontSize', fontsize, 'LineWidth', linewidth);
    bardat = [sm.pctTimeInRuns*100 baselinead.simpleMetrics.pctTimeInRuns*100];
    hhh = bar(1, bardat(1), 0.8, 'FaceColor', [1 1 1], 'EdgeColor', [0 0 0], 'LineWidth', linewidth); hold on
    hhh(2) = bar(2,  bardat(2), 0.8, 'FaceColor', [0 0 1], 'EdgeColor', [0 0 0], 'LineWidth', linewidth); hold off
    set(fracrunsaxes, 'FontSize', fontsize, 'FontName', font, 'YLim', [0 100], 'Box', 'off','XTick',[], 'LineWidth', linewidth);
    xl = get(fracrunsaxes, 'Xlim'); yl = get(fracrunsaxes, 'YLim');
    t1 = text (mean(xl), max(yl), 'time in runs', 'FontName', font, 'FontSize', fontsize, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
    p = get(t1, 'Extent');
    text (mean(xl), p(2), '(%)', 'FontName', font, 'FontSize', fontsize, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top');
else
    pos = [1-dw-7*w/8 rtpos(2)-h-spacing w*.75 h];
    rmpos = pos;
    fracrunsaxes = axes('Position', pos, 'Box', 'Off', 'FontName', font, 'FontSize', fontsize, 'LineWidth', linewidth);
    bardat = sm.pctTimeInRuns*100;
    hhh = bar(1, bardat, 0.8*2/3, 'FaceColor', [1 1 1], 'EdgeColor', [0 0 0], 'LineWidth', linewidth);
    set(fracrunsaxes, 'FontSize', fontsize, 'FontName', font, 'YLim', [0 100], 'Box', 'off','XTick',[], 'LineWidth', linewidth);
    xl = get(fracrunsaxes, 'Xlim'); yl = get(fracrunsaxes, 'YLim');
    t1 = text (mean(xl), max(yl), 'time in runs', 'FontName', font, 'FontSize', fontsize, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
    p = get(t1, 'Extent');
    text (mean(xl), p(2), '(%)', 'FontName', font, 'FontSize', fontsize, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top');
end

if (hasbaseline)
    pos = [dw rmpos(2)-h-spacing w h];
    lengthaxes = axes('Position', pos, 'Box', 'Off', 'FontName', font, 'FontSize', fontsize, 'LineWidth', linewidth);
    errorbar(lengthaxes,1, 10*sm.meanSpineLength, 10*sm.stdSpineLength, 'ko', 'LineWidth', linewidth); hold (lengthaxes, 'on');
    errorbar(lengthaxes,2, 10*baselinead.simpleMetrics.meanSpineLength, 10*baselinead.simpleMetrics.stdSpineLength, 'b*', 'LineWidth', linewidth); hold (lengthaxes, 'off')
    
    set(lengthaxes, 'FontSize', fontsize, 'FontName', font, 'YLim', [0 lengthylim], 'Box', 'off','XTick',2, 'XTickLabel', 'base', 'LineWidth', linewidth);
    xl = get(lengthaxes, 'Xlim'); yl = get(lengthaxes, 'YLim');
    t1 = text (mean(xl), max(yl), 'larval length', 'FontName', font, 'FontSize', fontsize, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
    p = get(t1, 'Extent');
    text (mean(xl), p(2), '(mm)', 'FontName', font, 'FontSize', fontsize, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top');
else
    pos = [dw+w/4 rmpos(2)-h-spacing w/2 h];
    lengthaxes = axes('Position', pos, 'Box', 'Off', 'FontName', font, 'FontSize', fontsize, 'LineWidth', linewidth);
    errorbar(lengthaxes,1, 10*sm.meanSpineLength, 10*sm.stdSpineLength, 'ko', 'LineWidth', linewidth); hold (lengthaxes, 'on');

    set(lengthaxes, 'FontSize', fontsize, 'FontName', font, 'YLim', [0 lengthylim], 'Box', 'off','XTick',2, 'XTickLabel', 'base', 'LineWidth', linewidth);
    xl = get(lengthaxes, 'Xlim'); yl = get(lengthaxes, 'YLim');
    t1 = text (mean(xl), max(yl), 'larval length', 'FontName', font, 'FontSize', fontsize, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
    p = get(t1, 'Extent');
    text (mean(xl), p(2), '(mm)', 'FontName', font, 'FontSize', fontsize, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top');
end

if (hasbaseline)
    pos = [0.5-w/2 rmpos(2)-h-spacing w h];
    fracvalaxes = axes('Position', pos, 'Box', 'Off', 'FontName', font, 'FontSize', fontsize, 'LineWidth', linewidth);
    bardat = [sm.fractionOfHSMarkedValid*100 baselinead.simpleMetrics.fractionOfHSMarkedValid*100];
    hhh = bar(1, bardat(1), 0.8, 'FaceColor', [1 1 1], 'EdgeColor', [0 0 0], 'LineWidth', linewidth); hold on
    hhh(2) = bar(2,  bardat(2), 0.8, 'FaceColor', [0 0 1], 'EdgeColor', [0 0 0], 'LineWidth', linewidth); hold off
    set(fracvalaxes, 'FontSize', fontsize, 'FontName', font, 'YLim', [0 100], 'Box', 'off','XTick',[], 'LineWidth', linewidth);
    xl = get(fracvalaxes, 'Xlim'); yl = get(fracvalaxes, 'YLim');
    t1 = text (mean(xl), max(yl), 'HS Validity (%)', 'FontName', font, 'FontSize', fontsize, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
 
else
    pos = [0.5-3*w/8 rmpos(2)-h-spacing w*.75 h];
    fracvalaxes = axes('Position', pos, 'Box', 'Off', 'FontName', font, 'FontSize', fontsize, 'LineWidth', linewidth);
    bardat = sm.fractionOfHSMarkedValid*100;
    hhh = bar(1, bardat(1), 0.8*2/3, 'FaceColor', [1 1 1], 'EdgeColor', [0 0 0], 'LineWidth', linewidth);
   
    set(fracvalaxes, 'FontSize', fontsize, 'FontName', font, 'YLim', [0 100], 'Box', 'off','XTick',[], 'LineWidth', linewidth);
    xl = get(fracvalaxes, 'Xlim'); yl = get(fracvalaxes, 'YLim');
    t1 = text (mean(xl), max(yl), 'HS Validity (%)', 'FontName', font, 'FontSize', fontsize, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
   
end

if (hasbaseline)
    pos = [1-dw-w rmpos(2)-h-spacing w h];
    nhsaxes = axes('Position', pos, 'Box', 'Off', 'FontName', font, 'FontSize', fontsize, 'LineWidth', linewidth);
    bardat = [sm.meanHSPerReo baselinead.simpleMetrics.meanHSPerReo];
    hhh = bar(1, bardat(1), 0.8, 'FaceColor', [1 1 1], 'EdgeColor', [0 0 0], 'LineWidth', linewidth); hold on
    hhh(2) = bar(2,  bardat(2), 0.8, 'FaceColor', [0 0 1], 'EdgeColor', [0 0 0], 'LineWidth', linewidth); hold off
    set(nhsaxes, 'FontSize', fontsize, 'FontName', font, 'YLim', [0 numhsylim], 'Box', 'off','XTick',[], 'LineWidth', linewidth);
    xl = get(nhsaxes, 'Xlim'); yl = get(nhsaxes, 'YLim');
    t1 = text (mean(xl), max(yl), 'Mean HS/Turn', 'FontName', font, 'FontSize', fontsize, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
else
    pos = [1-dw-7*w/8 rmpos(2)-h-spacing 3*w/4 h];
    nhsaxes = axes('Position', pos, 'Box', 'Off', 'FontName', font, 'FontSize', fontsize, 'LineWidth', linewidth);
    bardat = [sm.meanHSPerReo];
    hhh = bar(1, bardat(1), 0.8*2/3, 'FaceColor', [1 1 1], 'EdgeColor', [0 0 0], 'LineWidth', linewidth); 
    set(nhsaxes, 'FontSize', fontsize, 'FontName', font, 'YLim', [0 numhsylim], 'Box', 'off','XTick',[], 'LineWidth', linewidth);
    xl = get(nhsaxes, 'Xlim'); yl = get(nhsaxes, 'YLim');
    t1 = text (mean(xl), max(yl), 'Mean HS/Turn', 'FontName', font, 'FontSize', fontsize, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
end
