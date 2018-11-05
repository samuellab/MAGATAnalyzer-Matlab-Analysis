function axnew = cloneaxes( axold )
%function axnew = cloneaxes( axold )
%   creates a new axes with all the properties of the old axes
%{
props = set(axold);
fn = fieldnames(props);

axnew = axes();
inds = cellfun(@(x)~isempty(x), vals);
fn = fn(inds);
vals = vals(inds);
%}
axnew = axes('Position',get(axold,'Position'),...
           'XAxisLocation','top',...
           'YAxisLocation','right');
fn = {'XLim', 'YLim'};
invalid = {'Children','Title','XLabel','YLabel','ZLabel'};
[fn, ~] = setdiff(fn, invalid);
vals = get(axold, fn);

for j = 1:length(fn)
    %fn{j}
    %vals{j}
    set(axnew, fn{j}, vals{j});
end
%linkaxes([axold,axnew]);

