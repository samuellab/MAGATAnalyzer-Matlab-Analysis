function track = fromJava(jTr, ptType, trInd, loadImageByIndex, loadContour, camcalinfo, minpts)
%Used in ExFromJava to load a Track object from a TrackExtrationJava.Track
% C:\Users\Natalie\Documents\GitHub\Matlab-Track-Analysis\@Track\fromFile.m

% set flags/parameters
if (~exist ('loadImageByIndex', 'var'))
    loadImageByIndex = [];
end
if (~exist ('loadContour', 'var') || isempty (loadContour))
    loadContour = true;
end
if (~exist ('camcalinfo', 'var'))   
    camcalinfo = [];
end
existsAndDefault('minpts', 0);

% Check point type, create track, and set some so's if camcalinfo exists
if (isa(ptType, 'LarvaTrackPoint'))
    track = LarvaTrack();
    if (~isempty(camcalinfo)) %real points instead of camera points
        track.so.stop_speed_cut = 0.01;
        track.so.start_speed_cut = 0.015;
        track.so.curv_cut = 50;
    end
elseif (isa(ptType, 'MaggotTrackPoint'))
    track = MaggotTrack();
    if (~isempty(camcalinfo)) %real points instead of camera points
        track.so.stop_speed_cut = 0.01;
        track.so.start_speed_cut = 0.015;
        track.so.curv_cut = 50;
    end
else
    track = WormTrack();
end

%Set locinFile? (trackID is *indicative* of locInFile...)
track.locInFile = jTr.getTrackID();

% Determine nPts
track.npts = jTr.getNumPoints;

% If the track is going to get clipped from the experiment, don't bother with the rest
if (track.npts >= minpts) 
    % ~~ Ensure that loadImageByIndex.length = nPts
    if (length (loadImageByIndex) ~= track.npts)
        if (length(loadImageByIndex) == 1)
            loadImageByIndex = repmat(loadImageByIndex, [1 track.npts]);
        else
            if (~isempty(loadImageByIndex))
                disp(['length of loadImageByIndex = ' num2str(length(loadImageByIndex)) ' but track.npts = ' num2str(track.npts)]);
            end
            %By default, don't load any images
            loadImageByIndex = true([1 track.npts]);%%TEMP
        end
    end
    
    % ~~ Loop through nPts pts & add points
    pts(track.npts) = ptType;
    for i = 0:(track.npts-1)
        pts(i+1) = ptType.fromJava(jTr.getPoint(int32(i)), loadImageByIndex(i+1), loadContour, camcalinfo);
    end

    % ~~ Set points
    track.pt = pts;
    
    % ~~ Set startFrame&endFrame
    track.startFrame = jTr.getStart.getFrameNum;
    track.endFrame = jTr.getEnd.getFrameNum;
end


end