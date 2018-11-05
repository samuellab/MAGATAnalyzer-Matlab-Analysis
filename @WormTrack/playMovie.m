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
%   Axes, subplot(2,2,1:4) -- if passed in, we use these axes.  
%           if Axes has length 1, we only draw the movie part, not the data
%           fields
%   LabelTurns = true -- if false, we don't mark sharp turns on the graph
%   highlightInds = []; -- specify (interped) inds to draw over movie
%   highlightLineType = 'r.'; 

ptbuffer = 200;
delayTime = 0.1;
cc = track.expt.camcalinfo;
axisSize = max(size(track.pt(1).imData));
if (axisSize <= 0)
    axisSize = 20;
end
if (~isempty(cc))
    axisSize = axisSize * cc.realUnitsPerPixel;
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
Axes = [];
highlightInds = [];
highlightLineType = 'r.';
LabelTurns = true;
varargin = assignApplicable(varargin);
if (isempty(Axes))
    for j = 1:4
        Axes(j) = subplot(2,2,j);
    end
end

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
if (length(Axes) >= 4)
    datafields(track, sstart:send, Axes(2:4));
end
handles = [];
for j = inds
    ts1 = tic();
    %subplot(2,2,1); hold off; cla
    cla(Axes(1));
    pt(j).drawTrackImage(cc,'fid', fid, 'Axes', Axes(1), varargin{:}); 
    
    %pause
    
    hold (Axes(1), 'on');
    sstart = sind(j) - ptbuffer;
    send = sind(j) + ptbuffer;
    if (sstart < 1)
        sstart = 1;
    end
    if (send > length(sloc))
        send = length(sloc);
    end
    
    plot (Axes(1), sloc(1,sstart:send), sloc(2,sstart:send), 'b.-');
    plot (Axes(1), sloc(1,sind(j)), sloc(2,sind(j)), 'bo', 'MarkerSize', 5);
    if (~isempty(highlightInds))
         plot (Axes(1), sloc(1,highlightInds), sloc(2,highlightInds),highlightLineType);
    end
   % pause
    if (LabelTurns)
        rloc = sloc(:,revinds(revinds > sstart & revinds < send));
        oloc = sloc(:,omegainds(omegainds > sstart & omegainds < send));
        bloc = sloc(:,blipinds(blipinds > sstart & blipinds < send));
        rloc = rloc(:,abs(rloc(1,:) - loc(1,j)) < axisSize/2 & abs(rloc(2,:) - loc(2,j)) < axisSize/2);
        oloc = oloc(:,abs(oloc(1,:) - loc(1,j)) < axisSize/2 & abs(oloc(2,:) - loc(2,j)) < axisSize/2);
        bloc = bloc(:,abs(bloc(1,:) - loc(1,j)) < axisSize/2 & abs(bloc(2,:) - loc(2,j)) < axisSize/2);
    
        axes(Axes(1));
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
    end
    axis (Axes(1), [loc(1,j) + [-axisSize/2 axisSize/2], loc(2,j) + [-axisSize/2 axisSize/2]]);
    axis (Axes(1), 'equal'); 
    axis (Axes(1), [loc(1,j) + [-axisSize/2 axisSize/2], loc(2,j) + [-axisSize/2 axisSize/2]]);
    hold (Axes(1), 'off');
    set(Axes(1), 'XTick', [], 'YTick', []);
    if (~isempty(track.run))
        t = [];
        if (track.isrun(sind(j)))
            t = 'run';
        end

        if (~isempty(track.reorientation) && any([track.reorientation.inds] == sind(j)))
            t = 'reorientation';
        end
        title (Axes(1), t);
    end
    
   % pause
    if (length(Axes) >= 4)
        handles = updateCenter(handles, track, sind(j), sstart, send, Axes(2:4));
    end
    
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

function datafields(track, inds, Axes)
    fields = {'deltatheta','scovRatio','speed'};
    track.calculateDerivedQuantity(fields);
    mult = [rad2deg(1), 1, 1];
    for k = 1:3
        %subplot(2,2,k+1); hold off
        hold (Axes(k), 'off');
        plot (Axes(k), track.dq.eti(inds), mult(k)*track.dq.(fields{k})(inds), 'k.', 'MarkerSize', 5);
        hold (Axes(k), 'on');
     
    
        if (~isempty(track.run))
            runinds = find([track.run.endInd] > inds(1) & [track.run.startInd] < inds(end));
        else
            runinds = [];
        end
        for j = runinds;
            plot (Axes(k), track.dq.eti(track.run(j).inds), mult(k)*track.dq.(fields{k})(track.run(j).inds), 'm-','LineWidth',2);
        end
        title(Axes(k), fields{k});
    end

    spfields = {{'dthetaHiThresh', 'dthetaLoThresh','straightThetaThresh'},{}, {}};
    spcolors = {{'r-','g-','m-'},{}, {}};
    spmirror = [1 0 0];
    for k = 1:3
        %subplot(2,2,k+1);
        x = track.dq.eti(inds);
        for j = 1:length(spfields{k})
            f = spfields{k}{j};
            c = spcolors{k}{j};
            y = repmat (mult(k)*track.so.(f), size(x));
            plot (Axes(k), x,y,c);
            if (spmirror(k));
                plot (Axes(k), x,-y,c);
            end
        end
    end
end

function handles = updateCenter(handles, track, cind, start, stop, Axes)
    if ~isempty(handles)
        for j = 1:length(handles)
            delete(handles(j));
        end
    end
    fields = {'deltatheta','scovRatio','speed'};
    mult = [rad2deg(1), 1, 1];
    for k = 1:3
       % subplot(2,2,k+1); hold on
        %ih = ishold;
        handles(k) = plot (Axes(k), track.dq.eti(cind), mult(k)*track.dq.(fields{k})(cind), 'c.','MarkerSize',25);
        xlim(Axes(k), [min(track.dq.eti(start)) max(track.dq.eti(stop))]);
        %if (~ih)
        %    hold off
        %end
        if (k == 2)
            ylim(Axes(k), [1 4]);
        end
    end
end   




