function h = polarBarPlot (thetaInDeg, r, c, varargin)
%function h = polarBarPlot (thetaInDeg, r, c, varargin)
%
%c = color
%varargin: 'locOfZero', angle in degrees where 0 should be located
%
%anything you want passed to patch

existsAndDefault('c', 'b');
locOfZero = 0;
varargin = assignApplicable(varargin);

[tid, I] = sort(thetaInDeg + locOfZero);
r = r(I);
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
 x1 = r.*cosd(tid);
 y1 = r.*sind(tid);
 x2 = r.*cosd(tid([2:end 1]));
 y2 = r.*sind(tid([2:end 1]));
 
 x0 = zeros(size(x1));
 y0 = x0;

% xx = reshape([x0(:) x1(:) x2(:) x0(:)]',[],1);
% yy = reshape([y0(:) y1(:) y2(:) y0(:)]',[],1);
% size(xx)
xx = [xx{:}];
yy = [yy{:}];

h = patch(xx, yy, c, varargin{:});


end

