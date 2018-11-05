function ah = summaryFigure (wnstruct, varargin)


fignum = 1;
axeslabels = true;
equizoom = false;
scaleFactor = 1; %change to peak intensity/255
shortSigmaTime = 0;
longSigmaTime = .5;
figtitle = '';
varargin = assignApplicable(varargin);
if (isempty(wnstruct.btd))
    ah = [];
    return;
end
[adim, po] = blank8x10Figure(fignum, 'nrows', 5, 'topmargin', 1/11, 'bottommargin', 0.5/11);

redcolor = 'r';
bluecolor = 'b';
magentacolor = 'm';
cyancolor = 'c';
blackcolor = 'k';
tta_tr = [-10 5];
hsa_tr = [-5 5];
tsta_tr = [-20 5];
vktr = [0 20];

acchscolor = [0 183 178]/255;
rejhscolor = [200 0 50]/255;
lgturncolor = [255 163 0]/255;
smturncolor = [155 0 255]/255;


textaxespos = [adim.w5/2 adim.h0+0.2*(1-adim.h0) 1-adim.w5 (1-adim.h0)*0.8];
textaxes = axes('Position', textaxespos, 'Box', 'Off', 'XLim', [0 1], 'YLim', [0 1], 'XTick', [], 'YTick', []);
axis(textaxes, 'off');
if (isempty(figtitle))
    if (isfield(wnstruct, 'fname'))
        figtitle = fileparts(wnstruct.fname{1});
    else
        figtitle = fileparts(wnstruct.btd(1).fname);
    end
end
figtitle(figtitle == '_') = '-';
if (isfield(wnstruct, 'es'))
    es = wnstruct.es;
else
    es = [wnstruct.btd.es];
end
numan = round(sum([es.numAnimals]));
numturn = sum([es.numturns]);
numexp = length(es);
antime = sum([es.animalTime])/3600;
msg = [num2str(numexp), ' experiments; ' num2str(numan) ' larvae; ' num2str(numturn) ' turns; ' num2str(antime, 2), ' hours of data'];

t1 = text(0,0, msg, 'FontName', po.font, 'FontSize', po.bigfontsize-2, 'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom', 'Parent', textaxes);
p = get(t1, 'Extent');
t2 = text(0,p(4), figtitle, 'FontName', po.font, 'FontSize', po.bigfontsize, 'FontWeight', 'bold', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom', 'Parent', textaxes, 'Interpreter', 'None');



pos = [adim.lx3, adim.h0 - adim.h, adim.w3, adim.h];
an = 1;
ah(an).name = ['tta'];
ah(an).pos = pos;
ah(an).axes = axes('position', ah(an).pos, po.axesopts{:});
yd = scaleFactor*wnstruct.turn_start;
ydl = lowpass1D(yd, longSigmaTime/median(diff(wnstruct.taxis)));
if (shortSigmaTime > 0)
    yd = lowpass1D(yd, shortSigmaTime/median(diff(wnstruct.taxis)));
end

yd = scaleFactor*wnstruct.turn_start;

sf = sum(wnstruct.turn_start(wnstruct.taxis < 0 & wnstruct.taxis > - wnstruct.kernelTime).^2)/ sum(wnstruct.convkernels.^2) * length(wnstruct.convkernels) / nnz(wnstruct.taxis < 0 & wnstruct.taxis > - wnstruct.kernelTime);
sf = scaleFactor*sqrt(sf);

try
    yd_f10 = lowpass1D(scaleFactor*wnstruct.firstten.turn_start, longSigmaTime/median(diff(wnstruct.taxis)));
    yd_s10 = lowpass1D(scaleFactor*wnstruct.secondten.turn_start, longSigmaTime/median(diff(wnstruct.taxis)));
    ah(an).handles = plot (wnstruct.taxis, yd, 'k-', wnstruct.taxis,ydl, 'r-', wnstruct.taxis,yd_f10, 'g--',wnstruct.taxis,yd_s10, 'b--',0:(-wnstruct.kernelDt):(-wnstruct.kernelTime), wnstruct.convkernels*sf, 'm-',po.plotOptions{:});
   % legend ('all', 'lp all', 'lp 0-10', 'lp 10-20');
catch me
    ah(an).handles = plot (wnstruct.taxis, yd, 'k-', wnstruct.taxis,ydl, 'r-',0:(-wnstruct.kernelDt):(-wnstruct.kernelTime), wnstruct.convkernels*sf, 'm-', po.plotOptions{:});
end
xlim(tta_tr);
xlabel ('time rel. to turn start (s)'); ylabel('<intensity change>');
hold on; plot([min(xlim) max(xlim)], [0 0], 'k:', [0 0], [min(ylim) max(ylim)], 'k:');

%title ('turn triggered average');

an = an+1;
pos(1) = adim.cx3;
ah(an).name = ['tsa'];
ah(an).pos = pos;
ah(an).axes = axes('position', ah(an).pos, po.axesopts{:});
ah(an).handles = plot (wnstruct.taxis, scaleFactor*lowpass1D(wnstruct.bigturn_start_onehs, longSigmaTime/median(diff(wnstruct.taxis))), 'k-',wnstruct.taxis, scaleFactor*lowpass1D(wnstruct.smallturn_start_onehs, longSigmaTime/median(diff(wnstruct.taxis))),'k-',  po.plotOptions{:});
set(ah(an).handles(1), 'Color', lgturncolor);
set(ah(an).handles(2), 'Color', smturncolor);
yl = get(gca, 'YLim');
if (abs(yl(1)) > abs(yl(2)))
    loc = 'SouthWest';
else
    loc = 'NorthWest';
end
set(legend('big', 'small', 'Location', loc), 'FontName', po.font, 'FontSize', po.fontsize);
xlim(tta_tr);
xlabel ('time rel. to turn start (s)'); 
%title ('turn triggered average');

an = an+1;
pos(1) = adim.rx3;
ah(an).name = 'hsa';
ah(an).pos = pos;
ah(an).axes = axes('position', ah(an).pos, po.axesopts{:});
if (shortSigmaTime > 0)
    yd1 = lowpass1D(scaleFactor*wnstruct.first_acchs_start, shortSigmaTime/median(diff(wnstruct.taxis)));
    yd2 = lowpass1D(scaleFactor*wnstruct.first_rejhs_start, shortSigmaTime/median(diff(wnstruct.taxis)));
else
    yd1 = scaleFactor*wnstruct.first_acchs_start;
    yd2 = scaleFactor*wnstruct.first_rejhs_start;
end
ah(an).handles = plot (wnstruct.taxis, yd1, 'k-',wnstruct.taxis, yd2,'k-',  po.plotOptions{:});
set(ah(an).handles(1), 'Color', acchscolor);
set(ah(an).handles(2), 'Color', rejhscolor);
hold on; plot([min(xlim) max(xlim)], [0 0], 'k:', [0 0], [min(ylim) max(ylim)], 'k:');
xlabel('time (s) rel to h.s. start');
xlim(hsa_tr);

an = an+1;
pos(1) = adim.lx3; pos(2) = pos(2) - adim.dh;
ah(an).name = 'ens_stat';
ah(an).pos = pos;
ah(an).axes = axes('position', ah(an).pos, po.axesopts{:});

% cv = wnstruct.btd.behaviorTriggeredDataMatrix('all', '', wnstruct.linname, 0);
% cvt = wnstruct.btd.behaviorTriggeredDataMatrix('turn', 'start', wnstruct.linname, 0);
% hx = linspace(-3,3,30);
% h1 = histc(cv, binEdgesFromCenters(hx));
% h1 = h1(1:end-1)/sum(h1)/median(diff(hx));
% h2 = histc(cvt, binEdgesFromCenters(hx));
% h2 = h2(1:end-1)/sum(h2)/median(diff(hx));
% 
% m1 = mean(cv);
% s1 = std(cv);
% 
% cvt = (cvt-m1)/s1;
% dkl = 0.5*(var(cvt) - log(var(cvt)) + mean(cvt).^2 - 1);
hx = wnstruct.hist_x;
h1 = wnstruct.hist_all_val;
h2 = wnstruct.hist_turn_val;
h1 = h1/sum(h1)/median(diff(hx));
h2 = h2/sum(h2)/median(diff(hx));

dkl = wnstruct.dkl_est_gauss;
ah(an).handles = plot (hx, h1, hx, h2, po.plotOptions{:});
xlim([-3 3]);
set(ah(an).axes, 'YTick', [0 .25 .5]);
text(-2.9, .49, {'DKL = ', num2str(dkl,3)}, 'Parent', ah(an).axes, 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top', 'fontName', po.font, 'fontSize', po.fontsize);
%0.8*max(get(ah(an).axes, 'YLim'))
an = an+1;
pos(1) = adim.cx3;
ah(an).name = ['trvsconv'];
ah(an).pos = pos;
ah(an).axes = axes('position', ah(an).pos, po.axesopts{:});
ah(an).handles = shadedErrorPlot (wnstruct.lx, 60*wnstruct.tr_vs_conv, 60*wnstruct.tr_vs_conv_eb, 60*wnstruct.tr_vs_conv_eb);
hold on;
if (~isfield(wnstruct, 'rate_rog'))
    wnstruct.rate_rog = wnstruct.ratefun{1}(wnstruct.lx);
end
ah(an).handles = [ah(an).handles plot(wnstruct.lx, 60*wnstruct.rate_rog, 'r-', po.plotOptions{:})];
% set(ah(an).axes, 'YColor', tacolor{m});
ylabel(ah(an).axes, 'turn rate (min-1)');
xlabel(ah(an).axes, 'conv. value');
    
    
acchsdelta = wnstruct.meanChangeAfterAcceptedHeadsweep - wnstruct.meanChangeAfterHeadsweep;
rejhsdelta = wnstruct.meanChangeAfterRejectedHeadsweep - wnstruct.meanChangeAfterHeadsweep;
  
acchsdelta_eb = sqrt(wnstruct.meanChangeAfterAcceptedHeadsweep_eb.^2 + wnstruct.meanChangeAfterHeadsweep_eb.^2);
rejhsdelta_eb = sqrt(wnstruct.meanChangeAfterRejectedHeadsweep_eb.^2 + wnstruct.meanChangeAfterHeadsweep_eb.^2);

bigturndelta = wnstruct.meanChangeBeforeBigTurn - wnstruct.meanChangeBeforeTurn;
smallturndelta = wnstruct.meanChangeBeforeSmallTurn - wnstruct.meanChangeBeforeTurn;

bigturndelta_eb = sqrt(wnstruct.meanChangeBeforeBigTurn_eb.^2+wnstruct.meanChangeBeforeTurn_eb.^2);
smallturndelta_eb = sqrt(wnstruct.meanChangeBeforeSmallTurn_eb.^2+wnstruct.meanChangeBeforeTurn_eb.^2);

%handles = barweb(barvalues, errors, width, groupnames, bw_title, bw_xlabel, bw_ylabel, bw_colormap, gridstatus, bw_legend, error_sides, legend_type)
 
an = an+1;
pos = ah(an-1).pos; pos(1) = adim.rx3;
ah(an).pos = pos;
ah(an).name = 'changes for turn and hs';
ah(an).axes = axes('position', ah(an).pos, po.axesopts{:});
bardat = scaleFactor*[acchsdelta rejhsdelta ;smallturndelta bigturndelta];
ebdat = scaleFactor*[acchsdelta_eb rejhsdelta_eb ;smallturndelta_eb bigturndelta_eb];
hh = barweb_marc(bardat, ebdat, [] , {'HS', 'Turn Size'});
set(hh.bars(:,1), 'FaceColor', redcolor, 'edgecolor', redcolor);
set(hh.bars(:,2), 'FaceColor', 'w', 'edgecolor', redcolor);


if (isfield(wnstruct, 'ShannonInfo'))    
    an = an+1;
    pos(1) = adim.lx3; pos(2) = pos(2) - adim.dh; pos(3) = adim.w3;
    ah(an).name = 'hs bit rate';
    ah(an).pos = pos;
    ah(an).axes = axes('position', ah(an).pos, po.axesopts{:});
    plot (ah(an).axes, wnstruct.ShannonInfo.allhs.tcent, wnstruct.ShannonInfo.allhs.bit_rate, 'k'); 
    ah(an).axes(2) = axes('Position',get(ah(an).axes,'Position'),...
               'XAxisLocation','top',...
               'YAxisLocation','right',...
               'Color','none',...
               'XColor','b','YColor','b');
    plot (ah(an).axes(2),  wnstruct.ShannonInfo.allhs.tstart_cum, wnstruct.ShannonInfo.allhs.cumbit); 
    set(ah(an).axes(2), po.axesopts{:}, 'Color', 'none','XAxisLocation','top','YAxisLocation','right', 'XColor','b','YColor','b', 'XTickLabel', {}); 
    set(ah(an).axes, 'XLim', [-5 5]);
    xlabel(ah(an).axes(1), 't(s) rel to hs');
    ylabel(ah(an).axes(1), 'bits/s');
    %ylabel(ah(an).axes(2), 'bits');

    an = an+1;
    pos(1) = adim.cx3*.95 + adim.lx3*0.05; 
    ah(an).name = 'turn bit rate';
    ah(an).pos = pos;
    ah(an).axes = axes('position', ah(an).pos, po.axesopts{:});
    plot (ah(an).axes, wnstruct.ShannonInfo.turnsize.tcent, wnstruct.ShannonInfo.turnsize.bit_rate, 'k'); 
    ah(an).axes(2) = axes('Position',get(ah(an).axes,'Position'),...
               'XAxisLocation','top',...
               'YAxisLocation','right',...
               'Color','none',...
               'XColor','b','YColor','b');
    plot (ah(an).axes(2),  wnstruct.ShannonInfo.turnsize.tstart_cum, wnstruct.ShannonInfo.turnsize.cumbit); 
    set(ah(an).axes(2), po.axesopts{:}, 'Color', 'none','XAxisLocation','top','YAxisLocation','right', 'XColor','b','YColor','b', 'XTickLabel', {}); 
    set(ah(an).axes, 'XLim', [-30 15]);
    xlabel(ah(an).axes(1), 't(s) rel to turn');
    %ylabel(ah(an).axes(1), 'bits/s');
    ylabel(ah(an).axes(2), 'bits');
    
    
    an = an+1;
    pos = ah(an-1).pos; pos(1) = adim.rx3;
    ah(an).pos = pos;
    ah(an).name = 'changes for turn and hs';
    ah(an).axes = axes('position', ah(an).pos, po.axesopts{:});
    ind1a = find(wnstruct.ShannonInfo.allhs.tend_cum < 0, 1, 'last');
    ind1b = find(wnstruct.ShannonInfo.allhs.tstart_cum >= 0, 1, 'first');
    ind2a = find(wnstruct.ShannonInfo.turnsize.tend_cum < 0, 1, 'last');
    ind2b = find(wnstruct.ShannonInfo.turnsize.tstart_cum >= 0, 1, 'first');
    
    %ind1 = find(2*wnstruct.ShannonInfo.allhs.tcent-wnstruct.ShannonInfo.allhs.tstart >= 0, 1, 'first');
    %ind2 = find(2*wnstruct.ShannonInfo.turnsize.tcent-wnstruct.ShannonInfo.turnsize.tstart >= 0, 1, 'first');
    bardat = [wnstruct.ShannonInfo.allhs.totalbits wnstruct.ShannonInfo.allhs.cumbit(ind1a) wnstruct.ShannonInfo.allhs.totalbits-wnstruct.ShannonInfo.allhs.cumbit(ind1b);...
        wnstruct.ShannonInfo.turnsize.totalbits wnstruct.ShannonInfo.turnsize.cumbit(ind2a) wnstruct.ShannonInfo.turnsize.totalbits-wnstruct.ShannonInfo.turnsize.cumbit(ind2b)];
    hh = bar(bardat);
    set(hh(1), 'FaceColor', redcolor, 'edgecolor', redcolor);
    set(hh(2), 'FaceColor', 'w', 'edgecolor', redcolor);
    set(hh(3), 'FaceColor', 'w', 'edgecolor', 'k');
    set(ah(an).axes, 'XTickLabel', {'HS', 'Turn Size'})
   % axis tight
    yl = get(ah(an).axes, 'YLim');
    ylim(ah(an).axes, [0 yl(2)*1.1]);
else
    disp ('rerun analyzeBTDDirectory_onesource to get information calculations');
end
    
an = an+1;
pos(1) = adim.lx4; pos(2) = pos(2) - adim.dh; pos(3) = adim.w4;
ah(an).name = 'speed wiener';
ah(an).pos = pos;
ah(an).axes = axes('position', ah(an).pos, po.axesopts{:});

ah(an).handles = plot (-wnstruct.taxis, wnstruct.wiener.speed.k1);
xlim([0 -min(wnstruct.taxis)]);
set(ah(an).axes, 'XDir', 'reverse');
xlabel('tau (s)'); title ('k1 - speed');

an = an+1;
pos(1) = adim.clx4; 
ah(an).name = 'body bend wiener';
ah(an).pos = pos;
ah(an).axes = axes('position', ah(an).pos, po.axesopts{:});

ah(an).handles = plot (-wnstruct.taxis, wnstruct.wiener.absspineTheta.k1);
xlim([0 -min(wnstruct.taxis)]);
set(ah(an).axes, 'XDir', 'reverse');
xlabel('tau (s)'); title ('k1 - body bend mag');

an = an+1;
pos(1) = adim.crx4; 
ah(an).name = 'delta body bend wiener';
ah(an).pos = pos;
ah(an).axes = axes('position', ah(an).pos, po.axesopts{:});

ah(an).handles = plot (-wnstruct.taxis, wnstruct.wiener.dabsspineTheta.k1);
xlim([0 -min(wnstruct.taxis)]);
set(ah(an).axes, 'XDir', 'reverse');
xlabel('tau (s)'); title ('k1 - delta body bend mag');

an = an+1;
pos(1) = adim.rx4; 
ah(an).name = 'path dtheta wiener';
ah(an).pos = pos;
ah(an).axes = axes('position', ah(an).pos, po.axesopts{:});

ah(an).handles = plot (-wnstruct.taxis, wnstruct.wiener.absdeltatheta.k1);
xlim([0 -min(wnstruct.taxis)]);
set(ah(an).axes, 'XDir', 'reverse');
xlabel('tau (s)'); title ('k1 - path dtheta');

if (isfield(wnstruct, 'fname'))
    names = wnstruct.fname;
else
    names = {wnstruct.btd.fname};
end
for j = 1:length(names)
    [~,names{j}] = fileparts(names{j});
    names{j} = names{j}(end-13:end);
end
 pos(1) = adim.lx4; pos(2) = pos(2) - adim.dh*1.05; pos(3) = adim.w4;
 xp = [adim.lx4 adim.clx4 adim.crx4 adim.rx4];
for j = 1:4
    an = an+1;
    pos(1) = xp(j);
    ah(an).name = 'tta by exp';
    ah(an).pos = pos;
    ah(an).axes = axes('position', ah(an).pos, po.axesopts{:});
    inds = find(mod((1:length(wnstruct.exp)) - j,4) == 0);
    if (isempty(inds))
        continue;
    end
    ah(an).handles = plot (wnstruct.taxis, lowpass1D([wnstruct.turn_start; [wnstruct.exp(inds).turn_start]'], longSigmaTime/median(diff(wnstruct.taxis))));
    set(ah(an).handles(1), 'color', 'k', 'linestyle', ':');
   
    title (names(inds), 'interpreter', 'none');
    axis tight;
     xlim(tta_tr);
end
yl = get([ah((an-3):an).axes], 'YLim');
set([ah((an-3):an).axes], 'YLim', [min([yl{:}]) max([yl{:}])]);


% an = an+1;
% pos(1) = adim.rx4; 
% ah(an).name = 'path curv wiener';
% ah(an).pos = pos;
% ah(an).axes = axes('position', ah(an).pos, po.axesopts{:});
% 
% ah(an).handles = plot (-wnstruct.taxis, wnstruct.wiener.abscurv.k1);
% xlim([0 -min(wnstruct.taxis)]);
% set(ah(an).axes, 'XDir', 'reverse');
% xlabel('tau (s)'); title ('k1 - path curv');
% 
