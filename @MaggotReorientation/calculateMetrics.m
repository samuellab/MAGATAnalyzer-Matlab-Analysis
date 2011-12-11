function calculateMetrics(reo, prevRun, nextRun)
% calcuates several metrics (e.g. previous direction) and updates head sweep metrics too
% function calculateMetrics(reo)
% outputs: none
% inputs:
%   REO < MaggotReorientation
%   PREVRUN < Run, run that ended immediately before this reorientation
%       started
%   NEXTRUN < Run, run that starts immediately after this reorientation
%       ends

if (~isempty(prevRun) && isa(prevRun, 'Run'))
    startInd = prevRun.endInd + 1;
else
    startInd = [];
end
if (~isempty(nextRun) && isa(nextRun, 'Run'))
    endInd = max(startInd,nextRun.startInd - 1);
else
    endInd = [];
end


reo.numHS = length(reo.headSwing);
if (false && ~isempty(reo.headSwing))
    reo.startInd = min([reo.headSwing.startInd]);
    reo.endInd = max([reo.headSwing.endInd]);
else
    reo.startInd = startInd;
    reo.endInd = endInd;
end
reo.inds = reo.startInd:reo.endInd;

if (~isempty(reo.prevRun))
    reo.prevDir = reo.prevRun.endTheta;
end
if (~isempty(reo.nextRun))
    reo.nextDir = reo.nextRun.startTheta;
end

for j = 1:length(reo.headSwing)
    reo.headSwing(j).prevDir = reo.prevDir;
    reo.headSwing(j).nextDir = reo.nextDir;
end
