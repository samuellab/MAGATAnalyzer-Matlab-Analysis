function calculateMetrics(hs)
    if (length(hs) > 1)
        for j = 1:length(hs)
            hs(j).calculateMetrics;
        end
        return;
    end

    hs.track.calculateDerivedQuantity({'sbodytheta', 'shead', 'stail', 'smid', 'dbodytheta'});
    btheta = hs.track.dq.sbodytheta(hs.startInd:hs.endInd);
    
    [~,I] = max(abs(btheta));
    last = find(sign(hs.track.dq.sbodytheta(hs.startInd:hs.endInd)) == sign(btheta(I)), 1, 'last');
    if (~isempty(last))
         hs.endInd = hs.startInd + last - 1;
    end
    
    last = find(~hs.track.isrun(hs.startInd:hs.endInd),1, 'last');
    if (~isempty(last))
        hs.endInd = min(hs.endInd, hs.startInd + last);
    end
    if (hs.startInd > hs.endInd)
        disp('start > end')
        hs.startInd
        hs.endInd
    end
    %back up the start of the headswing to the point at which we passed the
    %head sweep start cutoff
    
    %first = find(abs(hs.track.dq.sbodytheta(1:hs.startInd)) < abs(hs.track.dq.sbodytheta(hs.startInd)) & sign(hs.track.dq.dbodytheta(1:hs.startInd)) == sign(hs.track.dq.sbodytheta(hs.startInd)), 1, 'last');
    first = find((abs(hs.track.dq.sbodytheta(1:hs.startInd)) < hs.track.so.headswing_start & sign(hs.track.dq.dbodytheta(1:hs.startInd)) == sign(hs.track.dq.sbodytheta(hs.startInd))) | sign(hs.track.dq.sbodytheta(1:hs.startInd)) ~= sign(hs.track.dq.sbodytheta(hs.startInd)), 1, 'last');
    if (~isempty(first))
        first = first + 1;
        %if we changed the start index, back up the previous run ending to
        %before the starting index as well
        if (first < hs.startInd && ~isempty(hs.prevRun)) %only back up the head swing start, do not move it forward
            %make sure the first head swing doesn't make the previous run too
            %short 
            first = max(first, ceil(hs.prevRun.startInd + hs.track.so.minRunTime/hs.track.dr.interpTime));
            first = min(first, hs.startInd);
            ei = min(hs.prevRun.endInd, first - 1);
            if (ei ~= hs.prevRun.endInd)
                hs.prevRun.endInd = ei;
                hs.prevRun.calculateMetrics(false);
            end
            hs.startInd = first;
        end
    end
    if (hs.startInd > hs.endInd)
        disp('2 start > end')
        hs.startInd
        hs.endInd
    end
    hs.inds = hs.startInd:hs.endInd;
   % btheta = hs.track.dq.sspineTheta(hs.inds); %changed from sbodytheta 7/2 by marc
    btheta = hs.track.dq.spineTheta(hs.inds); %changed from sspineTheta 12/10/2013 by marc
    
    mh = hs.track.dq.shead(:,hs.inds) - hs.track.dq.smid(:,hs.inds);
    tm = hs.track.dq.smid(:,hs.inds) - hs.track.dq.stail(:,hs.inds);
    
    
    [~,I] = max(abs(btheta));
    hs.maxInd = hs.inds(I);
    hs.maxTheta = btheta(I);
    hs.sign = sign(hs.maxTheta);
    hs.headDir = atan2(mh(2,I), mh(1,I));
    hs.tailDir = atan2(tm(2,I), tm(1,I));
    hs.accepted = hs.track.isrun(hs.endInd);
    hs.valid = all(hs.track.getDerivedQuantity('ihtValid', false, 'inds', hs.inds));
    if (~isempty(hs.nextRun))
        hs.nextDir = hs.nextRun.startTheta;
    end
    if (~isempty(hs.prevRun))
        hs.prevDir = hs.prevRun.endTheta;
    end
    %{
    hs.prevRun = hs.track.run(find([hs.track.run.stop] < hs.endInd, 1, 'last'));
    hs.nextRun = hs.track.run(find([hs.track.run.start] < hs.startInd, 1, 'first'));
    if (~isempty(hs.prevRun))
        hs.prevDir = hs.prevRun.endTheta;
    end    
    hs.nextDir
    %}
end