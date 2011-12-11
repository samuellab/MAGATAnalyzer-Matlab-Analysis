function cycleTime = detectCycleTime (expt,varargin)
%function cycleTime = detectCycleTime (expt, varargin)
%uses the fact that in current implementations, there is a 1 second delay
%whenever the slide switches to determine the cycle time
%
%the cycle time is defined as the number of frames per half period;  in
%other words how many frames of light there are at one time
%
%at this point, no optional arguments

inds = find(diff(expt.elapsedTime) > 0.75 & diff (expt.elapsedTime) < 1.25);
expt.cycleTime = median(diff(inds));

if nargout > 0
    cycleTime = expt.cycleTime;
end


