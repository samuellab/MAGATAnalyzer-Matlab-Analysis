function track = segmentTrack (track, mso)
%function track = segmentTrack (track, mso)
%

debug = false;

if (~exist('mso', 'var') || isempty(mso))
    mso = track.so;
else
    track.so = mso;
end

cv = track.getDerivedQuantity('curv');
bt = track.getDerivedQuantity('sbodytheta');
vdp = track.getDerivedQuantity('vel_dp');
sp = track.getDerivedQuantity(mso.speed_field);

if (debug)
    tx = track.dq.eti;
    figure(1);
    plot (tx, abs(cv), tx, repmat(mso.curv_cut, size(tx))); ylim([0 mso.curv_cut]);
    title ('curv')
    figure(2);
    plot (tx, bt, tx, repmat(mso.theta_cut, size(tx)));
    title ('body theta');
    figure(3);
    plot (tx, sp, tx, repmat(mso.stop_speed_cut, size(tx)),'r-', tx, repmat(mso.start_speed_cut, size(tx)),'g-');
    title ('speed');
end

highcurv = (abs(cv) > mso.curv_cut);
head_swinging = (abs(bt) > mso.theta_cut);
speedlow = (sp < mso.stop_speed_cut);

%whenever the head swings wide or the path has high 
%curvature or the speed drops too low, we say any existing run
%ends
notarun = (highcurv | head_swinging | speedlow);
endarun = find(diff(notarun) >= 1) + 1;


%in order to begin a run, we need to (a) not be at a stop point
%(b) be moving fast enough and (c) have the head aligned with the direction
%of motion
speedhigh = (sp >= mso.start_speed_cut);
headaligned = (vdp >= mso.aligned_dp); 
isarun = (~notarun & speedhigh & headaligned);

startarun = find(diff(isarun) >= 1) + 1;

start = startarun;
stop = endarun;
si = 1;
k = 0;
%create a list of sequential starts and stops from valid start and stop
%points
while (~isempty(si) && ~isempty(startarun))
    k = k+1;
    start(k) = startarun(si);
    ei = find(endarun > start(k), 1, 'first');
    if (isempty(ei))
        stop(k) = length(track.dq.eti);
    else
        stop(k) = endarun(ei);
    end
    si = find(startarun > stop(k),1, 'first');
end
start = start(1:k);
stop = stop(1:k);

%remove runs that are too short
inds = find(track.dq.eti(stop) - track.dq.eti(start) >= mso.minRunTime);
start = start(inds);
stop = stop(inds);
%{
%clear old runs
if (~isempty(track.run))
    for j = 1:length(track.run)
        delete(track.run(j));
    end  
end
%}
run = repmat(Run(),1);
%record runs in track, and take some basic statistics
for k = 1:length(start)
    run(k) = Run(track,start(k),stop(k));
   %{
    if (k > 1)
        run(k).previousRun = run(k-1);
        run(k-1).nextRun = run(k);
    end
    %}
end
track.run = run;
track.isrun = false(size(track.dq.eti));
track.isrun([run.inds]) = true;
notrun = ~track.isrun;

%eliminate all headswings before the first run & after the last run
firstrunind = find(track.isrun, 1, 'first');
lastrunind = find(track.isrun, 1, 'last');
inrange = false(size(notrun));
inrange(firstrunind:lastrunind) = true;

%now that we have found runs, we look for head swings
%locate head swings;  head swings are anything where the head swings 
head_swinging = find (abs(bt) > mso.headswing_start & notrun & inrange);
not_head_swing = find((abs(bt) < mso.headswing_stop) | ([0 diff(sign(bt))] ~= 0) & inrange);



%create a list of sequential starts and stops from valid start and stop
%points
si = 1;
k = 0;
start = head_swinging;
stop = start;
while (~isempty(si) &&~isempty(head_swinging))
    k = k+1;
    start(k) = head_swinging(si);
    ei = find(not_head_swing > start(k), 1, 'first');
    if (isempty(ei))
        stop(k) = length(track.dq.eti);
    else
        stop(k) = not_head_swing(ei);
    end
    si = find(head_swinging > stop(k), 1, 'first');
end
start = start(1:k);
stop = stop(1:k);


inds = start;
j = 0;
%a headswing is only valid if it includes at least one point that is not
%part of a run
for k = 1:length(start)
    if (any(notrun(start(k):stop(k))))
        j = j + 1;
        inds(j) = k;
    end
end
inds = inds(1:j);

start = start(inds);
stop = stop(inds);

headSwing = repmat(HeadSwing(), 0);
for k = 1:length(start)
    headSwing(k) = HeadSwing(track, start(k), stop(k));
end
track.headSwing = headSwing;

%group headswings into reorientations
%a reorientation is the period between runs, whether or not that contains
%any headswings;
%a reorientation is a group of 1 or more headswings that fall between the
%same runs
nextrun = zeros(size(track.headSwing));
for j = 1:length(track.headSwing)
    ind = find([track.run.startInd] > track.headSwing(j).startInd, 1, 'first');
    if (isempty(ind))
        nextrun(j) = length(track.run) + 1;
    else
        nextrun(j) = ind;
    end
end
k = 0;
reorientation = repmat(MaggotReorientation(), [1 length(track.run)-1]);
for j = 1:length(reorientation)
    inds = find(nextrun == j+1);
    if (isempty(inds))
        reorientation(j) = MaggotReorientation(track, [], track.run(j), track.run(j+1)); 
    else
        reorientation(j) = MaggotReorientation(track, track.headSwing(inds));
    end
end
track.reorientation = reorientation;
