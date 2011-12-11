classdef ESetCleaner
    %rules for cleaning up sets of Experiments 
    %clean(esc, eset) will use these rules to discard tracks that likely
    %contain invalid data, according to the rules below
    
    properties
        
        minHTValid = 0.95; %remove all tracks where mean(ihtValid) < minHTValid
        minSpeed = 0; %remove all tracks where mean(speed) < minSpeed
        minDist = 0; %remove all tracks where max(displacement) < minDist 
        minPts = 0; %remove all tracks with less than minPts pts
        pruneRect = []; %to be implemented later
        trimRect = []; % to be implemented later
        rpmCut = 1000; %remove if the phase angle accumulated by the track (in revolutions) is greater than the elapsed time of the track (in min) times rpmCut (rec 1)
        minRevCut = 3; %do not remove tracks for excessive revolution if less than this many complete revolutions
        
        askFirst = true; %whether to show report and verify cleaning prior to removing any tracks from eset
        showFigsInReport = true; %whether to show graphs in report or just display text on command line
        
    end
    
    methods
        [num2clean, reportstring] = getReport(esc, eset); %displays how many tracks will be removed & statistics of those removals
        
        clean(esc, eset); %removes tracks that likely contain invalid data
    end
    
end

