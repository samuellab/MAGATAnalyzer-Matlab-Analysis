function embiggen (ax,fs)
if (nargin < 1)
   ax = gca;
end
existsAndDefault('fs', 16);

set(ax,'FontSize', fs);
elems = {'XLabel', 'YLabel', 'Title'};
for j = 1:length(elems)
   set (get(ax,elems{j}), 'FontSize', fs);
end