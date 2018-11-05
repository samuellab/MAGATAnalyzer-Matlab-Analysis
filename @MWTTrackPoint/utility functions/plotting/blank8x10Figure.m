function [axesdims, po] = blank8x10Figure (fignum, varargin)
%function [axesdims, po] = blank8x10Figure (fignum, varargin)


hspace = .6/11;
nrows = 4;
leftmargin = .5/8;
rightmargin = .5/8;
topmargin = 1/11;
bottommargin = 1/11;
varargin = assignApplicable(varargin);


allaxeswidth = 1 - leftmargin - rightmargin;
if (existsAndDefault('fignum', []))
    f = figure(fignum);
else
    f = figure();
end
clf(f);
axesdims.fignum = f;
ss = get(0, 'ScreenSize');

screenwidth = ss(3); screenheight = ss(4);
figratio = 8/10;
r = 0.8;
figheight = screenheight*r;
figwidth = figheight*figratio;
figpos = round([(screenwidth - figwidth)/2, (screenheight-figheight)/2, figwidth, figheight]);
    
set(f, 'Position', figpos);
get(f, 'Position');  %calling get f,position here prevents a weird bug where text is incorrectly spaced later

set(f, 'PaperType', 'usletter', 'PaperPosition', [0.25 0.25 8 10], 'PaperOrientation', 'portrait', 'color', 'w', 'inverthardcopy', 'off');


axesdims.h0 = 1-topmargin;
allaxesheight = axesdims.h0-bottommargin;

wspace2 = 0.08 * allaxeswidth;
wspace3 = 0.08 * allaxeswidth;
wspace4 = 0.08 * allaxeswidth;
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
    po.lineWidth = 1;
    po.font = 'Arial';
    po.fontsize = 7;
    po.bigfontsize = 14;
    po.color = 'k';
    po.axesopts = {'FontName', po.font, 'FontSize', po.fontsize, 'LineWidth', po.lineWidth/2, 'box', 'off'};
    po.plotOptions = {'LineWidth', po.lineWidth};
    po.labelOptions = {'Interpreter', 'Tex', 'FontSize', po.fontsize};

end