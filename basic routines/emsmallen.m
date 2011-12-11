function emsmallen (ax,varargin)
if (nargin < 1)
   ax = gca;
end
FontSize = 10;
FontName = 'Arial';
elems = {'XLabel', 'YLabel', 'Title'};
varargin = assignApplicable(varargin);

set(ax,'FontSize', FontSize, 'FontName', FontName);
for j = 1:length(elems)
   set (get(ax,elems{j}),'FontSize', FontSize, 'FontName', FontName, varargin{:});
end