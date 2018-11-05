function flagOmegaTurnsAndReversals (track, wormSegmentOptions)
% flags sharp turns in a worm track
% function flagOmegaTurnsAndReversals (track, wormSegmentOptions)
%
% called by flagReorientations which is called by SegmentTrack
% usually not called directly by end user
%
% outputs: none
% inputs: 
%   TRACK: member of the track class (note MaggotTrack segmentation
%       proceeds differently)
%   WORMSEGMENTOPTIONS: segmentation options, see WormSegmentOptions class
%       for details
wso = wormSegmentOptions;

dt = track.getDerivedQuantity ('deltatheta');
sp = track.getDerivedQuantity('speed');
sharpTurns =  abs(dt) > wso.dthetaHiThresh;
%{ 
%old
%using imclose to link close together sharp turns
sz = size(sharpTurns);
sz(sz > 1) = wso.joinSTpts;
if (wso.joinSTpts > 1)
    sharpTurns = imclose(sharpTurns, ones(sz));
end
%}

%expand from regions of high delta theta until hitting a low delta theta
%combined with a higher speed
while(1)
    sz = size(sharpTurns);
    sz(sz > 1) = 3;
    newst = imdilate(sharpTurns, ones(sz)) & (abs(dt) > wso.dthetaLoThresh | sp < wso.speedEndSTThresh);
    if (all (newst == sharpTurns))
        break;
    else
        sharpTurns = newst;
    end
end

start = find (diff([0 sharpTurns]) > 0);
stop = find (diff([sharpTurns 0]) < 0);



if (isempty(start))
    start = 1;
end
if (isempty(stop))
    stop = length(track.dq.deltatheta);
end
%should never happen, but safety first
if (stop(1) < start(1))
    stop = stop(2:end);
end
if (start(end) > stop(end))
    start = start(1:(end-1));
end

for j = 1:length(start)
   sharpTurn(j) = WormSharpTurn(track, start(j), stop(j));    
end

track.sharpTurn = sharpTurn;
