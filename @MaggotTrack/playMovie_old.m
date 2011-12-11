function playMovie(track, varargin)
% one of several routines to play back track movies
%@MaggotTrack
%playMovie(track, varargin)
%enter options as pairs, caps matter
%options, with defaults
%
%ptbuffer = 1000;
%delayTime = 0.05;
%axisSize (size of image or 120 * mean speed)
%inds = 1:length(track.pt);
%startLoc = []; > if startLoc & stopLoc are not empty, we run the movie
%between these two points
%stopLoc = []; >
%startTime = []; if startTime and stopTime are not empty, re run the movie
%between these two times
%stopTime = [];
%pass 'fid', [] to not load images from disk


ptbuffer = 200;
delayTime = 0.1;
set(0,'DefaultTextInterpreter', 'Latex');
if (isempty(track.expt.camcalinfo))
    mf = 1;
else
    mf = track.expt.camcalinfo.realUnitsPerPixel;
end
axisSize = mf*max(size(track.pt(1).imData));
if (axisSize <= 0)
    try 
        track.expt.openDataFile();
        pttemp = track.expt.reloadPoint(track.pt(1));
        axisSize = mf*max(size(pttemp.imData));
    catch
    end
end


if (axisSize <= 0)
    %axisSize = 50;
    axisSize = 8 * mean(track.getDerivedQuantity('speed'));
end
iinds = [];
inds = 1:length(track.pt);
startLoc = [];
stopLoc = [];
startTime = [];
stopTime = [];
if (~isempty(track.expt))
    track.expt.openDataFile;
    fid = track.expt.fid;
else
    fid = [];
end
varargin = assignApplicable(varargin);
if (~isempty(iinds))
    inds = track.getDerivedQuantity('mapinterpedtopts', false,iinds);
end
if (~isempty(startLoc) && ~isempty(stopLoc))
    [~,s] = track.nearestPoint (startLoc);
    [~,e] = track.nearestPoint (stopLoc);
    if (s > e)
        inds = e:s;
    else
        inds = s:e;
    end
end
if (~isempty(startTime) && ~isempty(stopTime))
    pt = [track.pt];
    inds = find([pt.et] >= startTime, 1, 'first'):find([pt.et] <= stopTime, 1, 'last');
end
pt = [track.pt];
loc = [pt.loc];
sloc = track.getDerivedQuantity('sloc');
sind = track.getDerivedQuantity('mapptstointerped');
track.calculateDerivedQuantity({'sbodytheta', 'speed', 'vel_dp', 'dsbodytheta', 'spheadperp', 'sspineTheta'});

pt = [track.pt];
sstart = sind(1) - ptbuffer;
send = sind(end) + ptbuffer;
if (sstart < 1)
    sstart = 1;
end
if (send > length(sloc))
    send = length(sloc);
end
datafields(track, sstart:send);
handles = [];
ccinfo = track.expt.camcalinfo;
for j = inds
    ts1 = tic();
    subplot(2,2,1); hold off; cla
    pt(j).drawTrackImage(ccinfo,'fid', fid, varargin{:}); hold on
    sstart = sind(j) - ptbuffer;
    send = sind(j) + ptbuffer;
    if (sstart < 1)
        sstart = 1;
    end
    if (send > length(sloc))
        send = length(sloc);
    end
    
    plot (sloc(1,sstart:send), sloc(2,sstart:send), 'b.-');
    plot (sloc(1,sind(j)), sloc(2,sind(j)), 'bo', 'MarkerSize', 5);
    axis ([loc(1,j) + [-axisSize/2 axisSize/2], loc(2,j) + [-axisSize/2 axisSize/2]]);
    axis equal; 
    axis ([loc(1,j) + [-axisSize/2 axisSize/2], loc(2,j) + [-axisSize/2 axisSize/2]]);
    hold off
    set(gca, 'XTick', [], 'YTick', []);
    if (~isempty(track.run))
        t = [];
        if (track.isrun(sind(j)))
            t = [t 'run '];
        end
        if (~isempty(track.headSwing) && any([track.headSwing.inds] == sind(j)))
            t = [t 'headsweep '];
            I = find([track.headSwing.startInd] <= sind(j) & [track.headSwing.endInd] >= sind(j));
            if (~isempty(I))
                if (track.headSwing(I).accepted)
                    t = [t 'accepted '];
                else
                    t = [t 'rejected '];
                end
            end
        end
        if (~isempty(track.reorientation) && any([track.reorientation.inds] == sind(j)))
            t = [t 'reorientation '];
        end
        title (t);
    end
    
    handles = updateCenter(handles, track, sind(j), sstart, send);
    
    timeleft = delayTime - toc(ts1);
    if (timeleft > 0)
        pause(timeleft);
    else
        pause(0.001);
    end
        
end
%{
drawtime
tracktime
otherplotstime
%}
end

function datafields(track, inds)
    thetafield = 'sspineTheta';
    if (strcmpi (false && track.so.method,'new'))
        fields = {thetafield,track.so.speed_field,'spheadperp'};
        ftitles = {'Bondy Bend Angle', track.so.speed_field, 'vhead perp'};
        mult = [rad2deg(1), 1, 1];
    else
        fields = {thetafield,track.so.speed_field,'vel_dp'};
        ftitles = {'Bondy Bend Angle', track.so.speed_field, 'Velocity Dot Product'};
        mult = [rad2deg(1), 1, 1];
    end
    for k = 1:3
        subplot(2,2,k+1); hold off
        plot (track.dq.eti(inds), mult(k)*track.dq.(fields{k})(inds), 'k.', 'MarkerSize', 5); hold on
        if (~isempty(track.headSwing))
            hsinds = find([track.headSwing.endInd] > inds(1) & [track.headSwing.startInd] < inds(end));
        else
            hsinds = [];
        end
        for j = hsinds;
            if (track.headSwing(j).accepted)
                c = 'g-';
            else
                c = 'r-';
            end
            plot (track.dq.eti(track.headSwing(j).inds), mult(k)*track.dq.(fields{k})(track.headSwing(j).inds), c, 'LineWidth', 2);
        end
    
        if (~isempty(track.run))
            runinds = find([track.run.endInd] > inds(1) & [track.run.startInd] < inds(end));
        else
            runinds = [];
        end
        for j = runinds;
            plot (track.dq.eti(track.run(j).inds), mult(k)*track.dq.(fields{k})(track.run(j).inds), 'm-','LineWidth',2);
        end
        title(ftitles{k});
    end
    if (false && strcmpi (track.so.method,'new'))
        spfields = {{'theta_cut', 'headswing_start', 'headswing_stop'},{'stop_speed_cut', 'start_speed_cut'}, {}};
        spcolors = {{'m-','g-','r-'},{'r-', 'g-'}, {'m-'}};
        spmirror = [1 0 1];
    else
        spfields = {{'theta_cut', 'headswing_start', 'headswing_stop'},{'stop_speed_cut', 'start_speed_cut'}, {'aligned_dp'}};
        spcolors = {{'m-','g-','r-'},{'r-', 'g-'}, {'m-'}};
        spmirror = [1 0 0];
    end
    for k = 1:3
        subplot(2,2,k+1);
        x = track.dq.eti(inds);
        for j = 1:length(spfields{k})
            f = spfields{k}{j};
            c = spcolors{k}{j};
            y = repmat (mult(k)*track.so.(f), size(x));
            plot (x,y,c);
            if (spmirror(k));
                plot (x,-y,c);
            end
        end
    end
end

function handles = updateCenter(handles, track, cind, start, stop)
    thetafield = 'sspineTheta';
    if ~isempty(handles)
        for j = 1:length(handles)
            delete(handles(j));
        end
    end
    if (false && strcmpi (track.so.method,'new'))
        fields = {thetafield,track.so.speed_field,'spheadperp'};
         mult = [rad2deg(1), 1, 1];
    else
        fields = {thetafield,track.so.speed_field,'vel_dp'};
         mult = [rad2deg(1), 1, 1];
    end
   
    for k = 1:3
        subplot(2,2,k+1); hold on
        ih = ishold;
        handles(k) = plot (track.dq.eti(cind), mult(k)*track.dq.(fields{k})(cind), 'c.','MarkerSize',25);
        xlim([min(track.dq.eti(start)) max(track.dq.eti(stop))]);
        if (~ih)
            hold off
        end
    end
end   


function markHS(track, ind)
    hsind = find([track.headSwing.startInd] < ind & [track.headSwing.endInd] > ind);
    if (~isempty(hsind))
        subplot(2,2,1); hold on
        if (track.headSwing(hsind).accepted);
            c = 'g-';
            t = 'accepted headswing';
        else
            c = 'r-';
            t = 'rejected headswing';
        end
        plot (track.dq.shead(1, track.headSwing(hsind).inds), track.dq.shead(2, track.headSwing(hsind).inds), c);
        title (t);
        subplot(2,2,4); hold off
        x = track.dq.speed(track.headSwing(hsind).inds);
        y = rad2deg(track.dq.sbodytheta(track.headSwing(hsind).inds));
        plot (x, y, c, 'LineWidth', 2); hold on
        plot (repmat(track.so.start_speed_cut,size(y)),y, 'g--');
        plot (x,rad2deg(repmat(track.so.headswing_start, size(x))), 'g--');
        plot (x,rad2deg(repmat(track.so.headswing_stop, size(x))), 'r--');
        plot (x,-rad2deg(repmat(track.so.headswing_start, size(x))), 'g--');
        plot (x,-rad2deg(repmat(track.so.headswing_stop, size(x))), 'r--');
        
        title(t);
        xlabel ('speed');
        ylabel('theta');
    else
        runind = find([track.run.startInd] < ind & [track.run.endInd] > ind);
         if (~isempty(runind))
            subplot(2,2,1); title ('running');
            subplot(2,2,4); hold off; title ('running');
            x = track.dq.speed(track.run(runind).inds);
            y = rad2deg(track.dq.sbodytheta(track.run(runind).inds));
            plot (x, y, 'b-', 'LineWidth', 2); hold on
            plot (repmat(track.so.start_speed_cut,size(y)),y, 'g--');
            plot (repmat(track.so.stop_speed_cut,size(y)),y, 'r--');
            plot (x,rad2deg(repmat(track.so.headswing_start, size(x))), 'g--');
            plot (x,rad2deg(repmat(track.so.headswing_stop, size(x))), 'r--');
            plot (x,rad2deg(repmat(track.so.theta_cut, size(x))), 'm--');
            plot (x,rad2deg(-repmat(track.so.headswing_start, size(x))), 'g--');
            plot (x,rad2deg(-repmat(track.so.headswing_stop, size(x))), 'r--');
            plot (x,rad2deg(-repmat(track.so.theta_cut, size(x))), 'm--');
            xlabel ('speed');
            ylabel('theta');
         end
    end
    
end



