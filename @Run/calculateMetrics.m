function calculateMetrics(run, moveEndPosition)
%run.calculateMetrics()
%based on the track and indices already set, it finds a bunch of stuff
%for easy reading later
    if (~exist('moveEndPosition', 'var'))
        moveEndPosition = false;
    end
    run.track.calculateDerivedQuantity ({'eti', 'sloc', 'theta', 'pathLength', 'acc', 'fastTailSpeed'});
    
    if (moveEndPosition)
        dsp = diff(run.track.dq.speed);
        %this is a change to better determine turn start time
        testinds = floor(run.endInd - run.track.dr.smoothTime/run.track.dr.interpTime):ceil(run.endInd+run.track.dr.smoothTime/run.track.dr.interpTime);
        testinds = testinds(testinds > run.startInd+run.track.so.minRunTime/(2*run.track.dr.interpTime) & testinds < length(run.track.dq.eti)); %avoid run getting too short -- shouldn't be a problem ever
        testinds = testinds(testinds < run.endInd | dsp(testinds) < 0); %if it's still slowing down, extend past observed end of the run
        [~,I] = min(run.track.dq.fastTailSpeed(testinds)); %maggot often backs its tail up at the beginning of a turn
        run.endInd = testinds(1) + I - 2;
        %end change
    end
    
    ptbuffer = floor(min(1,run.track.so.minRunTime/2)/run.track.dr.interpTime);
    
    
    run.inds = run.startInd:run.endInd;
    startTheta = mean(unwrap(run.track.dq.theta(run.startInd + (1:ptbuffer))));
    endTheta = mean(unwrap(run.track.dq.theta(run.endInd + 1 - (1:ptbuffer))));
    run.startTheta = mod(startTheta+pi,2*pi) - pi;
    run.endTheta = mod(endTheta+pi,2*pi) - pi;
    
    displacement = run.track.dq.sloc(:,run.endInd) - run.track.dq.sloc(:,run.startInd);
    run.meanTheta = atan2(displacement(2), displacement(1));
    run.euclidLength = sqrt(sum((displacement).^2)); 
    run.pathLength = run.track.dq.pathLength(run.endInd) - run.track.dq.pathLength(run.startInd);
    run.runTime = run.track.dq.eti(run.endInd) - run.track.dq.eti(run.startInd);
end

