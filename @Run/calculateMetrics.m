function calculateMetrics(run)
%run.calculateMetrics()
%based on the track and indices already set, it finds a bunch of stuff
%for easy reading later

    ptbuffer = min(1,run.track.so.minRunTime/2)/run.track.dr.interpTime;
    run.track.calculateDerivedQuantity ({'eti', 'sloc', 'theta', 'pathLength'});
    run.inds = run.startInd:run.endInd;
    startTheta = mean(unwrap(run.track.dq.theta(run.startInd + (1:ptbuffer))));
    endTheta = mean(unwrap(run.track.dq.theta(run.endInd - (1:ptbuffer))));
    run.startTheta = mod(startTheta+pi,2*pi) - pi;
    run.endTheta = mod(endTheta+pi,2*pi) - pi;
    
    displacement = run.track.dq.sloc(:,run.endInd) - run.track.dq.sloc(:,run.startInd);
    run.meanTheta = atan2(displacement(2), displacement(1));
    run.euclidLength = sqrt(sum((displacement).^2)); 
    run.pathLength = run.track.dq.pathLength(run.endInd) - run.track.dq.pathLength(run.startInd);
    run.runTime = run.track.dq.eti(run.endInd) - run.track.dq.eti(run.startInd);
end

