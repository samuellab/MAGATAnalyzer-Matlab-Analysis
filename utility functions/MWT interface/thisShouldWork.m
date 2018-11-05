javaaddpath('/jvm/Chore.jar')
jsa = javaArray('java.lang.String',22'); jsa(1) = java.lang.String('-p'); jsa(2) = java.lang.String('0.02652'); jsa(3) = java.lang.String('-M'); jsa(4) = java.lang.String('2'); jsa(5) = java.lang.String('-t'); jsa(6) = java.lang.String('15'); jsa(7) = java.lang.String('--shadowless'); jsa(8) = java.lang.String('-S'); jsa(9) = java.lang.String('--plugin'); jsa(10) = java.lang.String('Reoutline::exp'); jsa(11) = java.lang.String('--plugin'); jsa(12) = java.lang.String('Respine'); jsa(13) = java.lang.String('--plugin'); jsa(14) = java.lang.String('SpinesForward'); jsa(15) = java.lang.String('--plugin'); jsa(16) = java.lang.String('Flux::gate::+E,892px,1185px,16mm,16mm'); jsa(17) = java.lang.String('--map'); jsa(18) = java.lang.String('-o'); jsa(19) = java.lang.String('speed,bias'); jsa(20) = java.lang.String('/data/kerrr/multi_paper/taxtap_main/20120319_192311.zip'); jsa(21) = java.lang.String('--target'); jsa(22) = java.lang.String('/data/kerrr/multi_paper/test');
jsp = javaArray('CustomComputation',4); jsp(1) = Reoutline(); jsp(2) = Respine(); jsp(3) = SpinesForward(); jsp(4) = Flux();
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
