function [axesdims, po] = blankPPTFigure (fignum, varargin)

figwidthin = 6;
hspace = .6;
nrows = 2;
leftmargin = .5;
rightmargin = .5;
topmargin = .5;
bottommargin = .5;
figratio = 4/3;
varargin = assignApplicable(varargin);

figheightin = figwidthin/figratio;
hspace = hspace/figheightin;
leftmargin = leftmargin/figwidthin;
rightmargin = rightmargin/figwidthin;
topmargin = topmargin/figheightin;
bottommargin = bottommargin/figheightin;


allaxeswidth = 1 - leftmargin - rightmargin;
if (exist ('fignum', 'var'))
    f = figure(fignum);
else
    f = figure();
end
clf(f);
axesdims.fignum = f;
ss = get(0, 'ScreenSize');

screenwidth = ss(3); screenheight = ss(4);

r = 0.6;
figheight = screenheight*r;
figwidth = figheight*figratio;
figpos = round([(screenwidth - figwidth)/2, (screenheight-figheight)/2, figwidth, figheight]);
    
set(f, 'Position', figpos);
get(f, 'Position');  %calling get f,position here prevents a weird bug where text is incorrectly spaced later

%set(f, 'PaperSize', [6.5 6*figratio+0.5], 'PaperPosition', [0.25 0.25 6 6*figratio], 'PaperOrientation', 'portrait', 'color', 'w', 'inverthardcopy', 'off');
set(f, 'PaperType', 'USLetter', 'PaperPositionMode', 'auto', 'PaperOrientation', 'landscape', 'color', 'w', 'inverthardcopy', 'off');


axesdims.h0 = 1-topmargin;
allaxesheight = axesdims.h0-bottommargin;

wspace2 = 0.12 * allaxeswidth;
wspace3 = 0.06 * allaxeswidth;
wspace4 = 0.06 * allaxeswidth;
wspace5 = 0.05 * allaxeswidth;

axesdims.h = (allaxesheight - (nrows-1)*hspace)/nrows;
axesdims.dh = axesdims.h + hspace;

axesdims.w2 = min((allaxeswidth - wspace2)/2);%, h*1.61803399);
axesdims.w3 = (allaxeswidth - 2*wspace3)/3;
axesdims.w4 = (allaxeswidth - 3*wspace4)/4;
axesdims.w5 = (allaxeswidth - 4*wspace5)/5;

centerx = 0.5*(1+leftmargin -rightmargin);

axesdims.lx2 = centerx - wspace2/2 - axesdims.w2;
axesdims.rx2 = centerx + (wspace2)/2;

axesdims.lx3 = centerx - wspace3 - 3*axesdims.w3/2;
axesdims.cx3 = centerx - axesdims.w3/2;
axesdims.rx3 = centerx + axesdims.w3/2 + wspace3;


axesdims.lx4 = centerx - 3*wspace4/2 - 2*axesdims.w4;
axesdims.clx4 = centerx - wspace4/2 - axesdims.w4;
axesdims.crx4 = centerx + (wspace4)/2;
axesdims.rx4 = centerx + 3*(wspace4)/2+axesdims.w4;

axesdims.lx5 = centerx - 2*wspace5 - 2.5*axesdims.w5;
axesdims.clx5 = centerx - wspace5 - 1.5*axesdims.w5;
axesdims.cx5 = centerx - axesdims.w5/2;
axesdims.crx5 = centerx + axesdims.w5/2 + (wspace5);
axesdims.rx5 = centerx + 1.5*axesdims.w5 + 2*(wspace5);

if nargout > 1
    po.lineWidth = 1.5;
    po.font = 'Arial';
    po.fontsize = 12;
    po.bigfontsize = 16;
    po.color = 'k';
    po.axesopts = {'FontName', po.font, 'FontSize', po.fontsize, 'LineWidth', po.lineWidth/2, 'box', 'off', 'Color', 'white'};
    po.plotOptions = {'LineWidth', po.lineWidth};
    po.labelOptions = {'Interpreter', 'Tex', 'FontSize', po.fontsize};

end