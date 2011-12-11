function detectCollisions(expt, maxdist)
%marks the iscollision field of each track if it comes too close to another
%track
%function detectCollisions(expt, maxdist)
%ouputs: none
%inputs:
%EXPT: a member of the experiment class
%MAXDIST: if two tracks are within this distance of each other AT THE SAME
%         POINT IN TIME, they are marked as in collision at that time

debug = false;

%initialize all collision states to false
for j = 1:length(expt.track) 
    expt.track(j).iscollision = false(size(expt.track(j).getDerivedQuantity('eti')));
end

%iterate through all pairs of tracks
for j = 1:length(expt.track)
    x = expt.track(j).getDerivedQuantity('sloc');
    t = expt.track(j).getDerivedQuantity('eti');
    for k = (j+1):length(expt.track)
        %find the position of maggot 2 at times corresponding to maggot 1
        x2 = expt.track(k).fieldAtTime('sloc', t);
        %find any where they are closer than maxdist
        d = sum((x - x2).^2, 1);
        c = (d < maxdist^2);
        if any(c)
            %mark points of collision as true in both tracks
            expt.track(j).iscollision = expt.track(j).iscollision | c;
            inds = expt.track(k).indsAtTime(t(c));
            if (~isempty(inds))
                expt.track(k).iscollision(inds) = true;
            end
            if (debug)
                clf;
                expt.track(j).plotPath('sloc', 'k-', 'highlightinds', c); hold on;
                expt.track(k).plotPath('sloc', 'b-', 'highlightinds', inds, 'highlightlinetype', 'g*');
                pause
            end
        end
    end
end
