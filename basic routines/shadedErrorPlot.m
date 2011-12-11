function h = shadedErrorPlot (xdata, ydata, eup, edown, c, varargin)
%function h = shadedErrorPlot (xdata, ydata, eup, edown, c, varargin) or
%c is the color;  
%h(1:j) is the patch; h(j+1:2*j) is the line
%
%optional param/value pairs, lineOptions, faceOptions 
lineOptions = {};
faceOptions = {};
varargin = assignApplicable(varargin);

existsAndDefault('c', 'b');
existsAndDefault('leg', {});
existsAndDefault('edown', []);
if (iscell(c))
    cc = c;
else
    cc = {c};
end
if (iscell(xdata))
    xd = xdata;
else
    xd = {xdata};
end
if (iscell(ydata))
    yd = ydata;
else
    yd = {ydata};
end
if (iscell(eup))
    ep = eup;
else
    ep = {eup};
end
if (iscell(edown))
    ed = edown;
else
    if (~isempty(edown))
        ed = {edown};
    else
        ed = ep;
    end
end
xmax = max([xd{:}]);
xmin = min([xd{:}]);
ymax = max([yd{:}] + [ep{:}]);
ymin = min([yd{:}] - [ed{:}]);

cla
set(gca, 'NextPlot', 'replace');
for j = 1:length(xd)
    c = cc{min(j, length(cc))};
    xdata = xd{j};
    ydata = yd{j};
    eup = ep{j};
    edown = ed{j};
    

    if (ischar(c))
        c = char2rgb(c);
    end

%existsAndDefault('edown', eup);
    xdata = xdata(:)';
    ydata = ydata(:)';
    eup = eup(:)';
    edown = edown(:)';
    [xdata,I] = sort(xdata);
    ydata = ydata(I);
    eup = eup(I);
    edown = edown(I);


    h(j) = plot (xdata, ydata, 'Color', c, lineOptions{:}, varargin{:}); hold on;
end
hold off;
set(gca, 'XLim', [xmin xmax], 'YLim', [ymin ymax]);
if (~isempty(leg))
    legend(h, leg, legendOptions{:});
end
for j = 1:length(xd)
    c = cc{min(j, length(cc))};
    xdata = xd{j};
    ydata = yd{j};
    eup = ep{j};
    edown = ed{j};
    

    if (ischar(c))
        c = char2rgb(c);
    end

%existsAndDefault('edown', eup);
    xdata = xdata(:)';
    ydata = ydata(:)';
    eup = eup(:)';
    edown = edown(:)';
    [xdata,I] = sort(xdata);
    ydata = ydata(I);
    eup = eup(I);
    edown = edown(I);
    ccc = get(gca, 'Color');
    hh(j) = patch([xdata xdata(end:-1:1)], [ydata+eup ydata(end:-1:1)-edown(end:-1:1)], 0.75*ccc + 0.25*c, 'LineStyle', 'none', faceOptions{:}, varargin{:});
end
%set(hh)
set(gca, 'Children', [h hh]);
h = [hh h];
%{
ih =  ishold(get(h(1), 'Parent'));
hold (get(h(1), 'Parent'), 'on');
%plot (xdata, ydata+eup, 'r-', xdata(end:-1:1), ydata(end:-1:1)-edown(end:-1:1), 'g-');
%pause

if (~ih)
    hold (gca, 'off');
end
%}
end

