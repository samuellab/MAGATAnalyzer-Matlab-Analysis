function h = shadedErrorPlotPolar (thdata, rdata, eup, edown, c, varargin)
%function h = shadedErrorPlotPolar (thdata, rdata, eup, edown, c, varargin) or
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
if (iscell(thdata))
    td = thdata;
else
    td = {thdata};
end
if (iscell(rdata))
    rd = rdata;
else
    rd = {rdata};
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
% 
% xmax = max([xd{:}]);
% xmin = min([xd{:}]);
% ymax = max([yd{:}] + [ep{:}]);
% ymin = min([yd{:}] - [ed{:}]);

rmax = max(rd{:} + ep{:});
cla
set(gca, 'NextPlot', 'replace');
[caxunder,cax] = polarBackground(rmax);
hold(cax, 'on');
for j = 1:length(td)
    c = cc{min(j, length(cc))};
    tdata = td{j};
    rdata = rd{j};
    eup = ep{j};
    edown = ed{j};
    
    tdata = tdata(:)';
    rdata = rdata(:)';
    eup = eup(:)';
    edown = edown(:)';
    
    [tdata,I] = sort(tdata);
    rdata = rdata(I);
    eup = eup(I);
    edown = edown(I);
    
    tdata = tdata([1:end 1]);
    rdata = rdata([1:end 1]);
    eup = eup([1:end 1]);
    edown = edown([1:end 1]);
    
    
    if (ischar(c))
        c = char2rgb(c);
    end
    
    xx = rdata.*cos(tdata);
    yy = rdata.*sin(tdata);
    
  
    h(j) = plot (cax, xx, yy, 'Color', c, lineOptions{:}, varargin{:}); 
    xx = [rdata + eup, rdata(end:-1:1)-edown(end:-1:1)].*cos([tdata tdata(end:-1:1)]);
    yy = [rdata + eup, rdata(end:-1:1)-edown(end:-1:1)].*sin([tdata tdata(end:-1:1)]);

    ccc = get(caxunder, 'Color');
    hh(j) = patch(xx, yy, 0.75*ccc + 0.25*c, 'LineStyle', 'none', faceOptions{:}, varargin{:}, 'Parent', cax);
end
if (~isempty(leg))
    legend(h, leg, legendOptions{:});
end
%set(hh)
set(cax, 'Children', [h hh]);
h = [hh h];

end

