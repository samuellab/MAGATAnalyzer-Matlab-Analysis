function clickMovie(track, varargin)
%plays a movie repeatedly between two points determined by mouse click
%function clickMovie(track, varargin)
%
%plots the path or segmentation of the track, then waits for user to zoom
%press enter to continue
%then waits for the user to click on points.  Plays movie between the first
%point and last point clicked
%
%inputs:
%TRACK: a member of the Track class
%VARARGIN: anything that can be passed to Track/playMovie
f = gcf;
figure(f);
clf(f);
if (isempty(track.run))
    track.plotPath('loc', 'bd-', 'MarkerSize', 2, 'LineWidth', 0.5); hold on
    track.plotPath('sloc', 'rd-', 'MarkerSize', 2, 'LineWidth', 0.5); hold off
else
    track.plotSegmentation();
end
pause
[x,y] = getpts();

clf(f);
while true
    figure(f);
    track.playMovie('startLoc', [x(1);y(1)], 'stopLoc', [x(end);y(end)], varargin{:});
end