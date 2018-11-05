function plotContourC(cc, varargin)
%function plotContourC(cc, varargin)
%
rect = [-Inf Inf -Inf Inf];
varargin = assignApplicable(varargin);

ind = 1;
ih = true;
while (ind < size(cc,2))
    nextind = cc(2,ind) + ind;
    x = cc(1,(ind+1):nextind);
    y = cc(2,(ind+1):nextind);
    if (any(x >= rect(1) & x <= rect(2) & y >= rect(3) & y <= rect(4)))
        h = plot (x, y, varargin{:});
        ax = get(h, 'Parent');
        ih = ih && ishold(ax);
        hold(ax, 'on');
    end
    ind = nextind+1;
end

if (~ih)
    hold(ax, 'off');
end