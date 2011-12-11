function cleanTracks (expt, minFrames, minDist)
%deletes short or stubby tracks; consider ESetCleaner instead
%function cleanTracks (expt, minFrames, minDist)
%expt.cleanTracks(minFrames, minDist)
%
%deletes all tracks that have fewer than minFrames pts
%and all tracks that extend less than minDist in any direction
%
%ouputs: none
%inputs:
%EXPT: a member of the experiment class
%MINFRAMES: any track with less than this number of points is deleted
%MINDIST: any track that extends less than minDist is deleted

    for j = length(expt.track):-1:1
        d(j) = maxExcursion(expt.track(j));
    end
    badtracks = find((d < minDist) | ([expt.track.npts] < minFrames));
    length(badtracks)
    inds = setdiff (1:length(expt.track), badtracks);
    length(inds)
    delete (expt.track(badtracks));

    length(expt.track)
    sum(isvalid(expt.track))
    
    expt.track = expt.track(inds);
    length(expt.track)
    sum(isvalid(expt.track))
        
end %cleanTracks



function d = maxExcursion(track)
    %function d = maxExcursion(track)

    loc = [track.pt.loc];
    x = loc(1,:);
    y = loc(2,:);

    if (length(x) > 100)
        n = floor(length(x)/100);
        x = x(1:n:end);
        y = y(1:n:end);
    end

    x = repmat(x,fliplr(size(x)));
    y = repmat(y,fliplr(size(y)));
    ds = (x-x').^2 + (y-y').^2;
    d = sqrt(max(max(ds)));
end