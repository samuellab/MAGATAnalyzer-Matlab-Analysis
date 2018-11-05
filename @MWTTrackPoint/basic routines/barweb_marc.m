function handles = barweb_marc(barvalues, errors, width, groupnames, bw_title, bw_xlabel, bw_ylabel, bw_colormap, gridstatus, bw_legend, error_sides, legend_type, barx)

%
% Usage: handles = barweb(barvalues, errors, width, groupnames, bw_title, bw_xlabel, bw_ylabel, bw_colormap, gridstatus, bw_legend, error_sides, legend_type)
%
% Ex: handles = barweb(my_barvalues, my_errors, [], [], [], [], [], bone,
% [], bw_legend, 1, 'axis')
%
% barweb is the m-by-n matrix of barvalues to be plotted.
% barweb calls the MATLAB bar function and plots m groups of n bars using the width and bw_colormap parameters.
% If you want all the bars to be the same color, then set bw_colormap equal to the RBG matrix value ie. (bw_colormap = [1 0 0] for all red bars)
% barweb then calls the MATLAB errorbar function to draw barvalues with error bars of length error.
% groupnames is an m-length cellstr vector of groupnames (i.e. groupnames = {'group 1'; 'group 2'}).  For no groupnames, enter [] or {}
% The errors matrix is of the same form of the barvalues matrix, namely m group of n errors.
% Gridstatus is either 'x','xy', 'y', or 'none' for no grid.
% No legend will be shown if the legend paramter is not provided
% 'error_sides = 2' plots +/- std while 'error_sides = 1' plots just + std
% legend_type = 'axis' produces the legend along the x-axis while legend_type = 'plot' produces the standard legend.  See figure for more details
%
% The following default values are used if parameters are left out or skipped by using [].
% width = 1 (0 < width < 1; widths greater than 1 will produce overlapping bars)
% groupnames = '1', '2', ... number_of_groups
% bw_title, bw_xlabel, bw_ylabel = []
% bw_color_map = jet
% gridstatus = 'none'
% bw_legend = []
% error_sides = 2;
% legend_type = 'plot';
%
% A list of handles are returned so that the user can change the properties of the plot
% handles.ax: handle to current axis
% handles.bars: handle to bar plot
% handles.errors: a vector of handles to the error plots, with each handle corresponding to a column in the error matrix
% handles.legend: handle to legend
%
%
% See the MATLAB functions bar and errorbar for more information
%
% Author: Bolu Ajiboye
% Created: October 18, 2005 (ver 1.0)
% Updated: Dec 07, 2006 (ver 2.1)
% Updated: July 21, 2008 (ver 2.3)

% Get function arguments
% Get function arguments
if nargin < 2
	error('Must have at least the first two arguments:  barweb(barvalues, errors, width, groupnames, bw_title, bw_xlabel, bw_ylabel, bw_colormap, gridstatus, bw_legend, barwebtype)');
elseif nargin == 2
	width = 1;
	groupnames = 1:size(barvalues,1);
	bw_title = [];
	bw_xlabel = [];
	bw_ylabel = [];
	bw_colormap = jet;
	gridstatus = 'none';
	bw_legend = [];
	error_sides = 2;
	legend_type = 'plot';
    barx = 1:size(barvalues, 1);
elseif nargin == 3
	groupnames = 1:size(barvalues,1);
	bw_title = [];
	bw_xlabel = [];
	bw_ylabel = [];
	bw_colormap = jet;
	gridstatus = 'none';
	bw_legend = [];
	error_sides = 2;
	legend_type = 'plot';
    barx = 1:size(barvalues, 1);
elseif nargin == 4
	bw_title = [];
	bw_xlabel = [];
	bw_ylabel = [];
	bw_colormap = jet;
	gridstatus = 'none';
	bw_legend = [];
	error_sides = 2;
	legend_type = 'plot';
    barx = 1:size(barvalues, 1);
elseif nargin == 5
	bw_xlabel = [];
	bw_ylabel = [];
	bw_colormap = jet;
	gridstatus = 'none';
	bw_legend = [];
	error_sides = 2;
	legend_type = 'plot';
    barx = 1:size(barvalues, 1);
elseif nargin == 6
	bw_ylabel = [];
	bw_colormap = jet;
	gridstatus = 'none';
	bw_legend = [];
	error_sides = 2;
	legend_type = 'plot';
    barx = 1:size(barvalues, 1);
elseif nargin == 7
	bw_colormap = jet;
	gridstatus = 'none';
	bw_legend = [];
	error_sides = 2;
	legend_type = 'plot';
    barx = 1:size(barvalues, 1);
elseif nargin == 8
	gridstatus = 'none';
	bw_legend = [];
	error_sides = 2;
	legend_type = 'plot';
    barx = 1:size(barvalues, 1);
elseif nargin == 9
	bw_legend = [];
	error_sides = 2;
	legend_type = 'plot';
    barx = 1:size(barvalues, 1);
elseif nargin == 10
	error_sides = 2;
	legend_type = 'plot';
    barx = 1:size(barvalues, 1);
elseif nargin == 11
	legend_type = 'plot';
    barx = 1:size(barvalues, 1);
elseif nargin == 12
    barx = 1:size(barvalues, 1);
end
if (size(barvalues, 1) == 1)
    barvalues = [barvalues;barvalues];
    barvalues(2,:) = NaN;
    errors = [errors;errors];
    errors(2,:) = NaN;
    onegroup = true;
    barx = 1:2;
else
    onegroup = false;
end
%xdata = 1:size(barvalues, 1);
xdata = barx;

change_axis = 0;
ymax = 0;
ymin = 0;
if size(barvalues,1) ~= size(errors,1) || size(barvalues,2) ~= size(errors,2)
	error('barvalues and errors matrix must be of same dimension');
end

if (~isempty(xdata) && size(barvalues,1) ~= length(xdata))
    barvalues = barvalues';
    errors = errors';
end

numgroups = size(barvalues, 1); % number of groups
numbars = size(barvalues, 2); % number of bars in a group

if isempty(width)
    width = 1;
end

% Plot bars
for j = 1:numgroups
    data = NaN(size(barvalues)); %was zeros
    data(j,:) = barvalues(j,:);
    handles.bars(j,:) = bar(xdata, data, width, 'edgecolor','k', 'linewidth', get(gca,'LineWidth')); 	hold on
end
% 
% data = NaN(size(barvalues));
% handles.bars(j+1,:) =  bar(xdata, data, width, 'edgecolor','k', 'linewidth', 2); 	hold on

if ~isempty(bw_colormap)
    colormap(bw_colormap);
else
    colormap(jet);
end
if ~isempty(bw_legend) && ~strcmp(legend_type, 'axis')
    handles.legend = legend(bw_legend, 'location', 'best', 'fontsize',12);
    legend boxoff;
else
    handles.legend = [];
end

% Plot erros
for i = 1:numbars
    for j = 1:numgroups
        xx = handles.bars(j,i).XData; %changed 3/21/2018
        x(j) = mean(xx(:,j));
%         xx =get(get(handles.bars(j,i),'children'), 'xdata');
%         if (isempty(xx))
%             continue;
%         end
%         x(j) = mean(xx([1 3],j));
    end
    if (exist('x', 'var') && ~isempty(x))
        handles.errors(i) = errorbar(x, barvalues(:,i), errors(:,i), 'k', 'linestyle', 'none', 'LineWidth', get(gca, 'LineWidth'));
        ymax = max([ymax; barvalues(:,i)+errors(:,i)]);
        ymin = min([ymin; barvalues(:,i)-errors(:,i)]);
    end
    
end

if error_sides == 1
    set(gca,'children', flipud(get(gca,'children')));
end
if (ymin ~= ymax)
    ylim([1.1*ymin ymax*1.1]);
end
if (onegroup)
    xlim ([0.5 1.5]);
else 
%     numgroups
%     xdata(1)
%     xdata(2)
%     xdata(end)
%     xdata(end-1)
%     1.5*xdata(1)-xdata(2)
%     xdata(end)*1.5-xdata(end-1)
    xlim([1.5*xdata(1)-0.5*xdata(2) xdata(end)*1.5-0.5*xdata(end-1)]);
    %xlim([0.5 numgroups+0.5]);
end

if strcmp(legend_type, 'axis')
    yl = get(gca, 'YLim');
    
    yloc = min(yl) * 1.002 - max(yl)*0.002;
    for i = 1:numbars
        xdata = get(handles.errors(i),'xdata');
        if (onegroup)
            ul = 1;
        else 
            ul = length(xdata);
        end
        for j = 1:ul
            %    try
            handles.leg(i,j) = text(xdata(j),  yloc, bw_legend(i), 'Rotation', 60, 'fontsize', 12, 'HorizontalAlignment', 'right');
            %   catch me
            %      disp(me.getReport)
            % xdata
            % ymax
            % bw_legend
            % end
        end
    end
    set(gca,'xaxislocation','top', 'box', 'on');
end

if ~isempty(bw_title)
    title(bw_title, 'fontsize',14);
end
if ~isempty(bw_xlabel)
    xlabel(bw_xlabel, 'fontsize',14);
end
if ~isempty(bw_ylabel)
    ylabel(bw_ylabel, 'fontsize',14);
end
set(gca, 'xticklabel', groupnames, 'box', 'on', 'ticklength', [0 0], 'xtick',1:numgroups, 'xgrid','off','ygrid','off');
if (onegroup)
    set(gca, 'XTick', 1);
end

if (change_axis == 1)
    xx = cell2mat(get(handles.errors, 'XData')');
    
    set(gca, 'XTick',  xx(1,:));
end
if ~isempty(gridstatus) && any(gridstatus == 'x')
    set(gca,'xgrid','on');
end
if ~isempty(gridstatus) && any(gridstatus ==  'y')
    set(gca,'ygrid','on');
end

handles.baseline = plot (get(gca, 'XLim'), [0 0], 'k-', 'LineWidth', get(gca, 'LineWidth'), 'Color', get(gca, 'XColor'));

handles.ax = gca;

% bz = handles.bars(:);
% for j = 1:length(bz)
%     xd = get(bz(j), 'XData');
%     yd = get(bz(j), 'YData');
%     xd = xd(isfinite(yd));
%     yd = yd(isfinite(yd));
%     set(bz(j), 'XData', xd, 'YData', yd);
% end

hold off
