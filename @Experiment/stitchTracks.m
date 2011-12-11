function stitchTracks (expt,  maxFrameInterval, maxDist, varargin)
%merges tracks that may represent the same animal
%function stitchTracks (expt,  maxFrameInterval, maxDist, varargin)
%
%inputs: 
%EXPT: a member of the Experiment class
%MAXFRAMEINTERVAL: the maximum number of frames that may be missing between
%   the end of track1 and the start of track2 to merge them
%MAXDIST: the maximum distance between the end of track1 and the start
%   of track2 to merge them
%optional arguments
%'interactive', true/false (default true)
%   if interactive, prompts user to approve all stitching first
%'stitchcollisions', true/false (default false)
%   if true, stitch tracks even if they are marked as colliding at the end
%   point
%'detectcollisions', true/false (default true)
%   if stitchcollisions is false and detect collisions is true, we run 
%   collision detection with separation maxDist before stitching
%'rectSize', (default maxDist*3), size of region surrounding stitch point shown on
%            map of tracks in and out

    %sort into ascending order by start frame
    %{
    sf = zeros(size(expt.track));
    for j = 1:length(expt.track)
        sf(j) = expt.track(j).pt(1).ind;
    end
    %}
    
    interactive = true;
    stitchcollisions = false;
    rectSize = maxDist * 3;
    varargin = assignApplicable(varargin);
    t = [expt.track];
    [~,I] = sort([t.startFrame]);
    expt.track = expt.track(I);

    if (~stitchcollisions)
        expt.detectCollisions(maxDist);
    end
    
    %j = 1;
    sf = [expt.track.startFrame];
    %work backwards so we only have to iterate once
    ts = tic;
    alreadyMerged = [];
    for j = length(expt.track):-1:1
        if (~stitchcollisions && expt.track(j).iscollision(end))
            continue; %don't stitch tracks that end in collisions
        end
        
        lastValid = j + find(sf((j+1):end) > expt.track(j).endFrame + maxFrameInterval, 1, 'first');
        
        firstValid = j + find(sf((j+1):lastValid) > expt.track(j).endFrame, 1, 'first');
        
        for k = setdiff(firstValid:lastValid, alreadyMerged);
            
            if (expt.track(j).precedes(expt.track(k), maxFrameInterval, maxDist) && ...
                 (stitchcollisions || ~expt.track(k).iscollision(1)))  %don't stitch tracks that begin in collision
                
                if (~interactive || getApproval(expt, j, k, rectSize))                    
                    alreadyMerged = mergeTracks (expt, j, k, alreadyMerged);
                end
                break; %only merge one track
            end
            
        end
       
        if (mod(j,100) == 0) 
            disp ([num2str(j /(length(expt.track) - j) * toc(ts)) ' s remaining']);
        end
    end
    inds = setdiff(1:length(expt.track), alreadyMerged);
    disp(['merged ' num2str(length(alreadyMerged)) ' tracks']);
    delete(expt.track(alreadyMerged));
    expt.track = expt.track(inds);
    
    
end %stitchTracks

function alreadyMerged = mergeTracks (expt, j, k, alreadyMerged)
    %disp('merging');
    expt.track(j).merge(expt.track(k));
    alreadyMerged = [alreadyMerged k];
end

function approved = getApproval (expt, j, k, rectSize, ptbuffer)
    ts = tic;
    existsAndDefault('rectSize', 20);
    existsAndDefault('ptbuffer', 10);
    figure(1); clf(1); hold on
    t1 = expt.track(j);
    t2 = expt.track(k);
    pt1 = expt.track(j).pt(end);
    pt2 = expt.track(k).pt(1);
    indsExpression = ['find(abs(track.getDerivedQuantity(''eti'') - ' num2str(pt1.et) ') < 5)'];
    expt.track.plotPath('iloc', 'k-', 'LineWidth', 0.2, 'Color', [0.7 0.7 0.7]);
    expt.track.plotPath('iloc', 'k--', 'indsExpression', indsExpression, 'LineWidth', 2);
    
    expt.executeTrackFunction('plotPath', 'loc', 'k-', 'LineWidth', 0.2, 'Color', [0.7 0.7 0.7]);
    expt.track(j).plotPath('loc', 'r-','LineWidth', 2);
    expt.track(k).plotPath('loc', 'g-','LineWidth', 2);
    loc = pt1.loc;
    axis([loc(1)-rectSize, loc(1) + rectSize, loc(2)-rectSize, loc(2)+rectSize]); 
    axis equal
    axis([loc(1)-rectSize, loc(1) + rectSize, loc(2)-rectSize, loc(2)+rectSize]); 
    
    fd = pt2.ind - pt1.ind;
    pd = pt1.distance(pt2);
    
    title({[expt.fname ' tracks ' num2str(j) ' & ' num2str(k) ],['frame diff = ' num2str(fd) ' point dist = ' num2str(pd)]},...
        'Interpreter','none');
    
    
    expt.openDataFile;
    figure(2); clf(2);
    ind1 = max(1,length([expt.track(j).pt]) - 10);
   
    ind2 = min(length(expt.track(k).pt),10);

    pt1 = expt.track(j).pt(ind1:end);
    pt2 = expt.track(k).pt(1:ind2);
    
    %find all points that are in target rectangle in time range
    disprect = [loc(1) - rectSize, loc(1) + rectSize, loc(2) - rectSize, loc(2) + rectSize];
    tinds = expt.track.passesThroughBox(disprect,...
        [pt1(1).et pt2(end).et]);
    t = expt.track(tinds);
    pt = [t.pt];
    ploc = [pt.loc];
    pin = [pt.ind];
    pinds = find(abs(ploc(1,:) - loc(1)) < rectSize & abs(ploc(2,:) - loc(2)) < rectSize ...
        & pin >= expt.track(j).pt(end).ind - ptbuffer & pin <= expt.track(k).pt(1).ind + ptbuffer);
    pt = expt.reloadPoint(pt(pinds));
    ploc = ploc(:,pinds);
    pin = pin(pinds);
    
   % disp(['getApproval setup took ' num2str(toc(ts)) ' s']);
   %rectSize
    while(1)
        figure(2);
        %{
        for j = 1:length(pt1)
            pt1(j).drawTrackImage([],'fid', expt.fid);
            title ('end of track 1');
            pause(0.1);
        end
        for j = 1:length(pt2)
            pt2(j).drawTrackImage([],'fid', expt.fid);
            title ('start of track 2');
            pause(0.1);
        end
        %}
      %  min(pin)
      %  max(pin)
        for j = min(pin):max(pin)
            cla();
            for k = find(pin == j)
                pt(k).drawTrackImage();
                hold on;
            end    
            indsExpression = ['find(track.getDerivedQuantity(''eti'') <= ' num2str(expt.elapsedTime(max(pin)) + 10) ' & track.getDerivedQuantity(''eti'') >= '  num2str(expt.elapsedTime(min(pin))-10) ')' ];
            tactive = [t.endFrame] >= j;
            t1.plotPath('iloc', 'y-', 'indsExpression', indsExpression, 'LineWidth', 4);
            t2.plotPath('iloc', 'y-', 'indsExpression', indsExpression, 'LineWidth', 4);
            indsExpression = ['find(track.getDerivedQuantity(''eti'') <= ' num2str(expt.elapsedTime(j)) ' & track.getDerivedQuantity(''eti'') >= '  num2str(expt.elapsedTime(min(pin))) ')' ];
            
            if (any(tactive))
                t(tactive).plotPath('iloc', 'c-', 'indsExpression', indsExpression, 'LineWidth', 2);
            end
            if (any(~tactive))
                t(~tactive).plotPath('iloc', 'r-', 'indsExpression', indsExpression, 'LineWidth', 2);
            end
            axis (disprect);
%            title ([num2str(sum(pin == j)) ' points should be displayed']);
            pause(0.1);
        end
   
        key = input ('stitch tracks? : y/n/a(ll remaining without asking)', 's');
        if (isempty(key))
            continue;
        end
        switch(lower(key(1)))
            case 'y'
                approved = true;
                return;
            case 'n'
                approved = false;
                return;
            case 'a'
                approved = true;
                evalin('caller', 'interactive = false;');
                return;
        end
    end
end