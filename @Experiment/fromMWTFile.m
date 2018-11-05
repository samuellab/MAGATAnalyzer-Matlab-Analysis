function expt = fromMWTFile (inputname, camcalinfo, varargin)
%loads an experiment from a bin file
%function expt = fromMWTFile (inputname, camcalinfo, varargin)
%this is a static method of the Experiment class (Experiment.fromMWTFile)
%
%outputs: 
%EXPT, a member of the experiment class
%inputs:
%INPUTFNAME: directory containing something choreography can read
%CAMCALINFO: camera calibration struct (ask Marc); pass empty ([]) to ignore
%   default: []
%optional args - 'choreography_command', 'pluginArray' -- ask Marc
%interpTime, (interpolation time in seconds) -- override default selection
%of interp time

if (~exist ('camcalinfo', 'var'))   
    camcalinfo = [];
end

choreography_command = '-t 30 -s 0.1 -M 1 --shadowless --plugin Reoutline::exp --plugin Respine::0.23::tapered=0.28,1,2 --plugin SpinesForward::rebias --minimum-biased 3mm -S --nanless';
pluginArray = javaArray('CustomComputation',3); pluginArray(1) = Reoutline(); pluginArray(2) = Respine(); pluginArray(3) = SpinesForward(); 
interpTime = [];
varargin = assignApplicable(varargin);
jsa = commandLineToJavaStringArray(choreography_command);
jsa(length(jsa)+1) = java.lang.String(inputname);
try
    chore = javaMethod('doEverything','Choreography',jsa,pluginArray,true);
    dances = chore.extractReasonableDancers();
catch me
    disp ('choreography error:');
    disp (me.getReport());
    expt = repmat(Experiment(), 0);
    return;
end

ts1 = tic();
expt = Experiment();

expt.fname = inputname;
expt.camcalinfo = camcalinfo;
expt.elapsedTime = double(chore.extractTimes())';
if (isempty(interpTime))
    interpTime = percentile(diff(expt.elapsedTime), 0.1);
end
for j = 1:length(dances)
    track(j) = MWTTrack.fromDance(dances(j), camcalinfo);
    track(j).dr.interpTime = interpTime;
    if (~isempty(camcalinfo)) %real points instead of camera points
        track(j).so.stop_speed_cut = 0.01;
        track(j).so.start_speed_cut = 0.015;
        track(j).so.curv_cut = 50;
    end
end
expt.track = track;
[expt.track.expt] = deal(expt);
disp(['conversion to matlab struct took: ' num2str(toc(ts1),2) ' s']);