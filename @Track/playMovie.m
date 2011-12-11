function playMovie(track, varargin)
% plays a movie of the track with annotation
% function playMovie(track, varargin)
% outputs, none
% inputs:
%   TRACK < Track;
%   VARARGIN:
%   enter options as pairs, caps matter
%   options, with defaults
%
%   ptbuffer = 200; -- how many points on either side of the current point
%       to plot in annotations
%   delayTime = 0.05; -- interframe delay
%   axisSize (size of image or 50) -- size of axes containing image
%   inds = 1:length(track.pt); -- inds to play
%   iinds = []; interped inds;  if passed, we find inds =
%               gdq(mapinterpedtopts,iinds)
%   startLoc = [];
%   stopLoc = []; if startLoc & stopLoc are both not empty, we run the movie
%       between these two points
%   startTime = [];
%   stopTime = []; if startTime & stopTime are both not empty, we run the
%       movie between these times
%   locField = 'sloc'; what field to use to plot the path over the movie
%   image


ptbuffer = 200;
delayTime = 0.1;

axisSize = max(size(track.pt(1).imData));
if (axisSize <= 0)
    axisSize = 50;
end
iinds = [];
inds = 1:length(track.pt);
startLoc = [];
stopLoc = [];
startTime = [];
stopTime = [];
track.expt.openDataFile;
fid = track.expt.fid;
locField = 'sloc';
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
pt = [track.pt];
if (~isempty(startTime) && ~isempty(stopTime))
    s = find([pt.et] >= startTime, 1, 'first');
    e = find([pt.et] <= stopTime, 1, 'last');
    inds = s:e;
end
loc = [pt.loc];
sloc = double(track.getDerivedQuantity(locField));
sind = track.getDerivedQuantity('mapptstointerped');
track.calculateDerivedQuantity({'speed', 'deltatheta', 'ddtheta'});
sstart = sind(1) - ptbuffer;
send = sind(end) + ptbuffer;


if (sstart < 1)
    sstart = 1;
end
if (send > length(sloc))
    send = length(sloc);
end
st = [track.sharpTurn];
if (~isempty(st))
    st = st([st.endInd] > sstart & [st.startInd] < send);
    revinds = [st([st.typeCode] > 0).centralInd];
    omegainds = [st([st.typeCode] < 0).centralInd];
    blipinds = [st([st.typeCode] == 0).centralInd];
else
    revinds = [];
    omegainds = [];
    blipinds = [];
end
datafields(track, sstart:send);
handles = [];
for j = inds
    ts1 = tic();
    subplot(2,2,1); hold off; cla
    pt(j).drawTrackImage([],'fid', fid, varargin{:}); 
    
    %pause
    
    hold on
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
    
   % pause
    
    rloc = sloc(:,revinds(revinds > sstart & revinds < send));
    oloc = sloc(:,omegainds(omegainds > sstart & omegainds < send));
    bloc = sloc(:,blipinds(blipinds > sstart & blipinds < send));
    rloc = rloc(:,abs(rloc(1,:) - loc(1,j)) < axisSize/2 & abs(rloc(2,:) - loc(2,j)) < axisSize/2);
    oloc = oloc(:,abs(oloc(1,:) - loc(1,j)) < axisSize/2 & abs(oloc(2,:) - loc(2,j)) < axisSize/2);
    bloc = bloc(:,abs(bloc(1,:) - loc(1,j)) < axisSize/2 & abs(bloc(2,:) - loc(2,j)) < axisSize/2);
    
    
    if (~isempty(rloc))
        text(rloc(1,:), rloc(2,:), 'R', 'Color', 'r', 'HorizontalAlignment', 'Center');
    end
    if (~isempty(oloc))
        text(oloc(1,:), oloc(2,:), '\Omega', 'Color', 'm', 'HorizontalAlignment', 'Center');
    end
    if (~isempty(bloc))
        text(bloc(1,:), bloc(2,:), 'b', 'Color', 'g', 'HorizontalAlignment', 'Center');
    end
    %plot (sloc(1,ri), sloc(2,ri), 'r.', sloc(1,oi), sloc(2,oi), 'm.', sloc(1,bi), sloc(2,bi), 'g.',...
     %   'MarkerSize', 10);
    
    % pause
     
    axis ([loc(1,j) + [-axisSize/2 axisSize/2], loc(2,j) + [-axisSize/2 axisSize/2]]);
    axis equal; 
    axis ([loc(1,j) + [-axisSize/2 axisSize/2], loc(2,j) + [-axisSize/2 axisSize/2]]);
    hold off
    set(gca, 'XTick', [], 'YTick', []);
    if (~isempty(track.run))
        t = [];
        if (track.isrun(sind(j)))
            t = 'run';
        end

        if (~isempty(track.reorientation) && any([track.reorientation.inds] == sind(j)))
            t = 'reorientation';
        end
        title (t);
    end
    
   % pause
    
    handles = updateCenter(handles, track, sind(j), sstart, send);
    
    timeleft = delayTime - toc(ts1);
    if (timeleft > 0)
        pause(timeleft);
    else
        pause(0.01);
    end
        
   % pause
end
%{
drawtime
tracktime
otherplotstime
%}
end

function datafields(track, inds)
    fields = {'deltatheta','scovRatio','speed'};
    track.calculateDerivedQuantity(fields);
    mult = [rad2deg(1), 1, 1];
    for k = 1:3
        subplot(2,2,k+1); hold off
        plot (track.dq.eti(inds), mult(k)*track.dq.(fields{k})(inds), 'k.', 'MarkerSize', 5); hold on
     
    
        if (~isempty(track.run))
            runinds = find([track.run.endInd] > inds(1) & [track.run.startInd] < inds(end));
        else
            runinds = [];
        end
        for j = runinds;
            plot (track.dq.eti(track.run(j).inds), mult(k)*track.dq.(fields{k})(track.run(j).inds), 'm-','LineWidth',2);
        end
        title(fields{k});
    end

    spfields = {{'dthetaHiThresh', 'dthetaLoThresh','straightThetaThresh'},{}, {}};
    spcolors = {{'r-','g-','m-'},{}, {}};
    spmirror = [1 0 0];
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
    if ~isempty(handles)
        for j = 1:length(handles)
            delete(handles(j));
        end
    end
    fields = {'deltatheta','scovRatio','speed'};
    mult = [rad2deg(1), 1, 1];
    for k = 1:3
        subplot(2,2,k+1); hold on
        ih = ishold;
        handles(k) = plot (track.dq.eti(cind), mult(k)*track.dq.(fields{k})(cind), 'c.','MarkerSize',25);
        xlim([min(track.dq.eti(start)) max(track.dq.eti(stop))]);
        if (~ih)
            hold off
        end
        if (k == 2)
            ylim([1 4]);
        end
    end
end   




