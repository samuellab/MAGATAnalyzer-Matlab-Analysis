function h = figureForPrinting(h)
%function h = figureForPrinting(h)

if (~exist('h','var') || isempty(h))
    h = figure();
else
    figure(h);
end

clf(h);

ratio = 4/3;

set(h,'units','inches');
p = get(h,'position');
p(3) = 4;
p(4) = 4/ratio;
set(h,'position', p);
set(h, 'units','pixels');
set(h, 'PaperUnits', 'inches');
set(h,'PaperPosition', [1 1 4 4/ratio]);
set(h,'PaperOrientation', 'portrait');
