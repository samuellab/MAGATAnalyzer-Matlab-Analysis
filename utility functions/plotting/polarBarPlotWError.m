function [h,heb] = polarBarPlotWError (thetaInDeg, r, ru, rl, c, varargin)
%function [h,heb] = polarBarPlot (thetaInDeg, r, ru, rl, c, varargin)
%
%c = color
%
% ru = upper error, if all < 0, not used 
% rl = lower error; if empty, same as upper error, if all < 0, not used
%varargin: 'locOfZero', angle in degrees where 0 should be located
%           'curveEB', [false]/true, whether to curve the error bars to
%           follow the data

%anything you want passed to patch

existsAndDefault('c', 'b');
locOfZero = 0;
curveEB = false;
varargin = assignApplicable(varargin);

[tid, I] = sort(thetaInDeg + locOfZero);
r = r(I);
if (size(tid, 1) > 1)
    tid = tid';
end

existsAndDefault('rl', ru);

for j = 1:length(tid)
    dt = median(diff(tid))/6;
    ti = tid(j) + (-dt:dt);
    if (curveEB)
        ri = interp1([(tid - 360) tid (tid+360)], [r r r], ti);
    else
        ri = r(j);
    end
    if (isempty(rl))
        continue;
    end
    if (length(rl) >= j && rl(j) >= 0)
        r1 = r(j) - rl(j);
        lowerbarx{j} = (ri - rl(j)) .* cosd(ti); %#ok<*AGROW>
        lowerbary{j} = (ri - rl(j)) .* sind(ti);
        lowerstemx{j} = [r1 r(j)] * cosd(tid(j));
        lowerstemy{j} = [r1 r(j)] * sind(tid(j));
    else
        r1 = r(j);
        lowerbarx{j} = [];
        lowerbary{j} = [];
        lowerstemx{j} = [];
        lowerstemy{j} = [];
    end
    if (length(ru) >= j && ru(j) >= 0)
        r2 = r(j) + ru(j);
        upperbarx{j} = (ri + ru(j)) .* cosd(ti);
        upperbary{j} = (ri + ru(j)) .* sind(ti);
        upperstemx{j} = [r(j) r2] * cosd(tid(j));
        upperstemy{j} = [r(j) r2] * sind(tid(j));
    else
        r2 = r(j);
        upperbarx{j} = [];
        upperbary{j} = [];
        upperstemx{j} = [];
        upperstemy{j} = [];
    end
    
    
    
end

ti = min(tid):max(tid);
ri = interp1(tid, r, ti);
tid = tid - median(diff(tid))/2; %was diff (thetaInDegrees)
for j = 1:length(tid)
    
    if (j < length(tid))
        inds = (ti >= tid(j) & ti <= tid(j+1));
    else
        inds = [find(ti >= tid(j)) find(ti + 360 <= tid(j) + median(diff(tid))/2)];
    end
    xx{j} = [0 ri(inds).*cosd(ti(inds)) 0];
    yy{j} = [0 ri(inds).*sind(ti(inds)) 0];
end



% 
% 
% 
%  x1 = r.*cosd(tid);
%  y1 = r.*sind(tid);
%  x2 = r.*cosd(tid([2:end 1]));
%  y2 = r.*sind(tid([2:end 1]));
%  
%  x0 = zeros(size(x1));
%  y0 = x0;

% xx = reshape([x0(:) x1(:) x2(:) x0(:)]',[],1);
% yy = reshape([y0(:) y1(:) y2(:) y0(:)]',[],1);
% size(xx)
xx = [xx{:}];
yy = [yy{:}];

h = patch(xx, yy, c, varargin{:});
hold on;
if (isempty(rl) || (all(rl) < 0 && all(ru) < 0))
    heb = [];
    return;
end
for j = 1:length(upperstemx)
    heb(j,:) = [plot(upperstemx{j}, upperstemy{j}, upperbarx{j}, upperbary{j}, 'Color', c, 'LineWidth', 1);plot(lowerstemx{j}, lowerstemy{j}, lowerbarx{j}, lowerbary{j}, 'Color', 'w', 'LineWidth', 1)];
end


end

