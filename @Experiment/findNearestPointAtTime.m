function [pt, track, trackind, ptind] = findNearestPointAtTime(expt, loc, atTime, varargin)

    % Location/time 
    withinDist = 1; % cm within loc
    withinTime = 10; % seconds within time below
    
    varargin = assignApplicable(varargin);
        
    % Get the point
    if (~exist ('loc', 'var') || isempty (loc))
        [x,y] = getpts();   
        %use last point selected
        x = x(end);
        y = y(end);
    else
        x = loc(1,end);
        y = loc(2,end);
    end
    
    % Gather points from each track
    trackInds = [];
    pts = [];
    dist = [];
    for i=1:length(expt.track)
        [newPts, newDs] = gatherPointsInRange(expt.track(i),x,y, atTime, withinDist, withinTime);
        if ~isempty(newPts)
            pts = [pts newPts];
            dist = [dist newDs];
            newInds = i*ones(1,length(newPts));
            trackInds = [trackInds newInds];
        end
    end
    
    if ~isempty(dist)
        [~,sortInds] = sort(dist);
        pt = pts(sortInds(1));
        trackind = trackInds(sortInds(1));
        track = expt.track(trackind);
        ptind = pt.ind;
    else
        pt = [];
        trackind = [];
        track = [];
        ptind = [];
    end
    
    
end



function [pts, ds] = gatherPointsInRange(track, x,y, time, withinDist, withinTime)

    timeRange = [time-withinTime time+withinTime];
    if timeRange(1)<track.pt(1).et
        timeRange(1)=track.pt(1).et;
    end
    if timeRange(2)>track.pt(end).et
        timeRange(2)=track.pt(end).et;
    end
    
    % Check if the track has any proper time pts
    if (track.pt(end).et<timeRange(1) || track.pt(1).et>timeRange(2));
        pts = [];
        ds = [];
        return;
    end
    
    % Gather pts in time range
    tpp=(track.pt(end).et-track.pt(1).et)/length(track.pt);
    ptRange = floor((timeRange-track.pt(1).et)./tpp);
    if ptRange(1)<1
        ptRange(1)=1;
    end
    if ptRange(2)>length(track.pt)
        ptRange(2)=length(track.pt);
        
    end
    
    pts = track.pt(ptRange(1):ptRange(end));
    locs = [pts.loc];
    
    if ~isempty(pts)
        % Measure dist to those pts, return points within range
        dist = (locs(1,:)-x).^2 + (locs(2,:)-y).^2;

        pts = pts(dist<withinDist);
        ds = dist(dist<withinDist);
    else
        ds=[];
    end
end


