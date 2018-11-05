function [lineHandle, textHandle ] = labeledBar(ax, x, y, str, offset, c, lineOptions, textOptions )
%function [lineHandle, textHandle ] = labeledBar(ax, x, y, str, offset, c, lineOptions, textOptions )
%   creates a labeled line using coordinates from selected axis, using annotation

existsAndDefault('lineOptions', {});
existsAndDefault('textOptions', {});


[lx,ly] = dsxy2figxy_marc(ax, x, y);
lineHandle = annotation('line', lx, ly, 'Color', c, lineOptions{:});
switch (lower(offset))
    case 'right'
        tx = x(2); ty = mean(y);
        va = 'middle'; ha = 'left';
    case 'left'
        tx = x(1); ty = mean(y);
        va = 'middle'; ha = 'right';
   case 'above'
        tx = mean(x); ty = y(2);
        va = 'bottom'; ha = 'center';
    case 'below'
        tx = mean(x); ty = y(1);
        va = 'top'; ha = 'center';
    otherwise
        tx = mean(x); ty = mean(y); 
        va = 'middle'; ha = 'center';
end
textHandle = annotation('textbox', dsxy2figxy_marc(ax, [tx ty 0 0]), 'String', str, 'Color', 'w', 'LineStyle', 'none', 'EdgeColor', 'none', 'BackgroundColor', 'none', 'FitBoxToText', 'on', 'VerticalAlignment', va, 'HorizontalAlignment', ha, 'Margin', 2, textOptions{:});
end

% tp = textHandle.Position;
% lp = lineHandle.Position;
% switch (lower(offset))
%     case 'right'
%         lpm = lp(2) + 0.5 * lp(4);
%         tp(1) = lp(1) + lp(3);
%         tp(2) = lpm - tp(4)/2;
%         textHandle.Position = tp;
%     case 'left'
%         lpm = lp(2) + 0.5 * lp(4);
%         tp(1) = lp(1) - tp(3);
%         tp(2) = lpm - tp(4)/2;
%         textHandle.Position = tp;
%    case 'above'
%         lpm = lp(1) + 0.5 * lp(2);
%         tp(2) = lp(2) + lp(4);
%         tp(1) = lpm - tp(3)/2;
%         textHandle.Position = tp;
%     case 'below'
%         lpm = lp(1) + 0.5 * lp(2);
%         tp(2) = lp(2) - tp(4);
%         tp(1) = lpm - tp(3)/2;
%         textHandle.Position = tp;
%     otherwise
%         disp('unknown offset');
% end
%         
% 
% 
% end

