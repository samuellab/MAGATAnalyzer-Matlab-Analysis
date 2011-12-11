function addCycleTime (expt, varargin)
%expt.addCycleTime -- adds information about the cycle to the experiment
%and its tracks
%
%dq fields added are 
%lighton:  a boolean that is true if the light is on
%timeon: the time in seconds since the light last turned on;
%timeoff: the time in seconds since the light last turned off;

period = expt.cycleTime;
f = 1:length(expt.elapsedTime);

%first frame is off; second frame is on 
%period + 1 frame is on; period + 2 frame is off
lighton = mod(f-2, 2*period) < period;

oninds = find(diff(lighton) > 0) + 1;
offinds = find(diff(lighton) < 0) + 1;
timeon = expt.elapsedTime;
for j = length(oninds):-1:1
    ind1 = oninds(j);
    if (j < length(oninds))
        ind2 = oninds(j+1) - 1;
    else
        ind2 = length(timeon);
    end
    timeon(ind1:ind2) = timeon(ind1:ind2) - timeon(ind1);
end
timeon(1:oninds(1)) = NaN;

timeoff = expt.elapsedTime;
for j = length(offinds):-1:1
    ind1 = offinds(j);
    if (j < length(offinds))
        ind2 = offinds(j+1) - 1;
    else
        ind2 = length(timeoff);
    end
    timeoff(ind1:ind2) = timeoff(ind1:ind2) - timeoff(ind1);
end
timeoff(1:offinds(1)) = NaN;

gq = GlobalQuantity();
gq.xField = 'eti';
gq.derivationMethod = @GlobalQuantity.oneDinterpolation;
gq.xData = expt.elapsedTime;

gq.fieldname = 'lighton';
gq.yData = lighton;
expt.addGlobalQuantity(gq);

gq.fieldname = 'timeon';
gq.yData = timeon;
expt.addGlobalQuantity(gq);

gq.fieldname = 'timeoff';
gq.yData = timeoff;
expt.addGlobalQuantity(gq);

