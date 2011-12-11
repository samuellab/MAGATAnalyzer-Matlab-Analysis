function draw(run, varargin) 
%function draw(run, varargin) 
%
if (isempty(run.track))
    return;
end
Color = 'b-';
LineWidth = 1.5; 

varargin = assignApplicable(varargin);

ih = ishold;
plot (run.track.dq.sloc(1,run.inds), run.track.dq.sloc(2,run.inds), Color, 'LineWidth', LineWidth, varargin{:});
hold on;
t = rad2deg(run.track.dq.theta(run.inds(1)));
t = mod(t,360);
if (t < 45 || t >= 315)
    symbol = 'g>';
end
if (45 <=t && t < 135)
    symbol = 'g^';
end
if (135 <=t && t < 225)
    symbol = 'g<';
end
if (225 <= t)
    symbol = 'gv';
end

plot (run.track.dq.sloc(1,run.inds(1)), run.track.dq.sloc(2,run.inds(1)), symbol);
plot (run.track.dq.sloc(1,run.inds(end)), run.track.dq.sloc(2,run.inds(end)), 'rh');

if (~ih)
    hold off
end
