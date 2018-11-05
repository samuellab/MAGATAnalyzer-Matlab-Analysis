function calculateMetrics(reo)
%WormReorientation/
%function calculateMetrics(reo)
%
if (length(reo) > 1)
    for j = 1:length(reo)
        reo(j).calculateMetrics;
    end
    return;
end

if (isempty(reo) || isempty(reo.sharpTurn))
    return;
end

reo.startInd = reo.sharpTurn(1).startInd;
reo.endInd = reo.sharpTurn(end).endInd;
if (isempty(reo.prevRun))
    reo.thetaIn = reo.sharpTurn(1).thetaIn;
else 
    reo.thetaIn = reo.prevRun.endTheta;
end
if (isempty(reo.nextRun))
    reo.thetaOut = reo.sharpTurn(end).thetaOut;
else
    reo.thetaOut = reo.nextRun.startTheta;
end

reo.dTheta = diff(unwrap([reo.thetaIn;reo.thetaOut]));
  
reo.inds = reo.startInd:reo.endInd;

reo.numTurns = length(reo.sharpTurn);
reo.turnsequence = [reo.sharpTurn.typeCode];
