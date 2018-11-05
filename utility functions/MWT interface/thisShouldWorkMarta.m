javaaddpath('Chore.jar')
%cl = '-t 30 -s 0.1 -M 1 -p $PV --shadowless --plugin Reoutline::exp --plugin Respine::0.23::tapered=0.28,1,2 --plugin SpinesForward::rebias --minimum-biased 3mm -S --nanless';
cl = '-t 30 -s 0.1 -M 1 --shadowless --plugin Reoutline::exp --plugin Respine::0.23::tapered=0.28,1,2 --plugin SpinesForward::rebias --minimum-biased 3mm -S --nanless';

inputfname = '/Users/gershow/mwt worm data/for Marc/20130111_115604';
outputdir = '/Users/gershow/mwt worm data/for Marc/20130111_115604/test';
jsa = commandLineToJavaStringArray(cl);
jsa(length(jsa)+1) = java.lang.String(inputfname);
%jsa(length(jsa)+1) = java.lang.String('--target');
%jsa(length(jsa)+1) = java.lang.String(outputdir);

if (~exist(outputdir, 'dir'))
    mkdir(outputdir);
end

jsp = javaArray('CustomComputation',3); jsp(1) = Reoutline(); jsp(2) = Respine(); jsp(3) = SpinesForward(); 
%jsp(4) = Flux();
chore = javaMethod('doEverything','Choreography',jsa,jsp,true);
dances = chore.extractReasonableDancers();
length(dances)
xy = dances(1).extractNthSpinePoints(0);
plot(xy)
xy2 = dances(1).extractNthSpinePoints(10);
I = 1:(length(xy)/2);
plot(I,xy(I),I,xy2(I))



% Try also: extractCentroidPoints, extractOutlineAtFrame(n) -- counts up from zero
% Output arrays are x0 x1 x2 x3 .. xN y0 y1 y2 ... yN
