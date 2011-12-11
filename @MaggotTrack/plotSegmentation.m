function plotSegmentation (track, varargin)
%function plotSegmentation (track, varargin)
track.calculateDerivedQuantity({'sloc', 'stail','shead','theta'});
ih = ishold;
plot (track.dq.sloc(1,:), track.dq.sloc(2,:), 'm.-', 'LineWidth', 1); hold on
if (isempty(track.run))
    return;
end
for j = 1:length(track.run)
    track.run(j).draw
end

for j = 1:length(track.reorientation)
    track.reorientation(j).draw();
end
if (~ih)
    hold off
end

axis equal