function loadJavaMasks( track, jTr )
%LOADJAVAMASKS Summary of this function goes here
%   Detailed explanation goes here

    for i=1:length(track.pt)
        %This just adds a field which may or may not be accessible by other
        %functions, since it's not currently in the header
        track.pt(i).mask = jTr.getPoint(i-1).getMask().getIntArray()';%This may throw an error bc getPoint returns a point, not a mtp
    end
    
end

