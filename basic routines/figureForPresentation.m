function h = figureForPresentation(h)
%function h = figureForPrinting(h)

if (~exist('h','var') || isempty(h))
    h = figure();
else
    figure(h);
end

clf(h);

xdim = 1024;
ydim = 768;

set(h,'units','pixels');
p = get(h,'position');
p(3) = xdim;
p(4) = ydim;
set(h,'position', p);
set(h, 'units','pixels');
set(h, 'PaperUnits', 'inches');
set(h,'PaperPosition', [1 1 4 4/ratio]);
set(h,'PaperOrientation', 'portrait');
