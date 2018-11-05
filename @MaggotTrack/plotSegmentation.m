function plotSegmentation (track, varargin)
%function plotSegmentation (track, varargin)
TrackColor = 'm.-';
RunColor = 'b-';
ReoColors.accepted = 'g-';
ReoColors.rejected = 'r-';
SegColors = [];
varargin = assignApplicable(varargin);
if (isempty(SegColors))
    SegColors.TrackColor = TrackColor;
    SegColors.RunColor = RunColor;
    SegColors.ReoColors = ReoColors;
end

track.calculateDerivedQuantity({'sloc', 'stail','shead','theta'});
ih = ishold;
plot (track.dq.sloc(1,:), track.dq.sloc(2,:), SegColors.TrackColor, 'LineWidth', 1); hold on
if (isempty(track.run))
    return;
end
for j = 1:length(track.run)
    track.run(j).draw('Color', SegColors.RunColor);
end

for j = 1:length(track.reorientation)
    track.reorientation(j).draw('Colors', SegColors.ReoColors);
end
if (~ih)
    hold off
end

axis equal