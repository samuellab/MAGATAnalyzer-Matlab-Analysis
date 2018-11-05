function hh = rateFunVAdaptGraphs (rf, c, marker, fitType, varargin)
ax = gca;
varargin = assignApplicable(varargin);
ih = ishold(ax);

if (length(rf) > 1)    
    existsAndDefault('c', {[0 0 0]});
    existsAndDefault('fitType', 'linear');
    existsAndDefault('marker', 'o');
    for j = 1:length(rf)
        if (iscell(c) && length(c) >= j)
            cc = c{j};
        else
            cc = c;
        end
        if (iscell(fitType) && length(fitType) >= j)
            ff = fitType{j};
        else
            ff = fitType;
        end
         if (iscell(marker) && length(marker) >= j)
            mm = marker{j};
        else
            mm = marker;
         end
        if (j == 1)
            h = rateFunVAdaptGraphs (rf(j), cc, mm, ff, 'ax', ax, varargin{:});
            hold (ax, 'on');
        else
            h = [h;rateFunVAdaptGraphs(rf(j), cc, mm, ff, 'ax',ax, varargin{:})];
        end
        
    end
    if (nargout > 0)
        hh = h;
    end
    if (~ih)
        hold(ax, 'off');
    end
    return;
end
existsAndDefault('c', {[0 0 0]});
existsAndDefault('fitType', 'linear');
existsAndDefault('marker', 'o');
if (iscell(c))
    c = c{1};
end
if (iscell(fitType))
    fitType = fitType{1};
end
if (iscell(marker))
    marker = marker{1};
end


nclip = 0;
ax = gca;
rscale = 60;
varargin = assignApplicable(varargin);

h = errorbar(ax, rf.lx((1+nclip):(end-nclip)), rscale*rf.rate((1+nclip):(end-nclip)), rscale*rf.rate_eb((1+nclip):(end-nclip)), marker, 'Color', c); %#ok<*AGROW>
hold on;

switch (fitType)
    case 'both'
        hp = plot(rf.lx, rscale*exp(polyval(rf.rateFitLin, rf.lx)), 'k--', rf.lx, rscale*exp(polyval(rf.rateFitQuad, rf.lx)), 'k-');
    case 'quad'
        hp = plot(rf.lx, rscale*exp(polyval(rf.rateFitQuad, rf.lx)), 'k-');
    case 'none'
        hp = repmat(h(1),0);
    otherwise
        hp = plot(rf.lx, rscale*exp(polyval(rf.rateFitLin, rf.lx)), 'k--');
end
set(hp, 'Color', c);
h = [h;hp];
if (~ih)
    hold(ax, 'off');
end
if (nargout > 0)
    hh = h;
end