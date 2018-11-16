function playMovie(track, varargin)
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


ptbuffer = ceil(40/track.dr.interpTime);
delayTime = 0.05;
%set(0,'DefaultTextInterpreter', 'Latex');
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
vidObj = [];
aviobj = [];
avirect = [];
resizefig = true;
figcolor = 'w';
fontcolor = 'k';
if (~isempty(track.expt))
    track.expt.openDataFile;
    fid = track.expt.fid;
else
    fid = [];
end
if (isa(track.pt(1), 'MWTTrackPoint'))
    imOptions = {'pretty', true, 'drawHeadArrow', true, 'spineColor', 'y.','contourColor', 'm-', 'contourWidth', 1, 'scale', 1, 'drawContour', true, 'drawSpine', true};%'drawContour', true, 'scale', 1,'contourColor', 'r-', 'LineWidth', 1, 'mhWidth', 1};
else
    imOptions = {'pretty', true, 'drawHeadArrow', true, 'spineColor', 'y.','contourColor', 'k.', 'contourWidth', 1, 'scale', 1, 'drawContour', false, 'drawSpine', false};%'drawContour', true, 'scale', 1,'contourColor', 'r-', 'LineWidth', 1, 'mhWidth', 1};
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
sloc = double(track.getDerivedQuantity('sloc'));
sind = double(track.getDerivedQuantity('mapptstointerped'));
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

[imax, dataaxes] = setupFigure(gcf, resizefig);

u = get(imax, 'units');
set(imax, 'units', 'pixels');
imh = get(imax, 'position');
imh = imh(4);
fontsize = .04 * imh;
set(imax, 'units', u);
%set(gcf, 'DefaultAxesColor', 'w', 'DefaultTextColor', 'k');
set(gcf, 'Color', figcolor);
set(imax, 'Color', 'k');
set(dataaxes, 'XColor', fontcolor, 'YColor', fontcolor, 'FontSize', fontsize-2)
datafields(track, sstart:send, dataaxes, fontsize);


handles = [];
ccinfo = track.expt.camcalinfo;
eti = track.getDerivedQuantity('eti');
for j = inds
    ts1 = tic();
    hold (imax,'off'); 
    cla(imax); 
    set(imax, 'Color', 'k');
    pt(j).drawTrackImage(ccinfo,'fid', fid, 'Axes', imax, imOptions{:} ); hold (imax, 'on');
    sstart = sind(j) - ptbuffer;
    send = sind(j) + ptbuffer;
    if (sstart < 1)
        sstart = 1;
    end
    if (send > length(sloc))
        send = length(sloc);
    end
    
    plot (imax,sloc(1,sstart:send), sloc(2,sstart:send), 'w.-', 'MarkerSize', 5);
%    plot (imax, sloc(1, sstart:(1/track.dr.interpTime):send), sloc(2, sstart:(1/track.dr.interpTime):send), 'w.');
    plot (imax,sloc(1,sind(j)), sloc(2,sind(j)), 'wo', 'MarkerSize', 5);
    axis (imax,[loc(1,j) + [-axisSize/2 axisSize/2], loc(2,j) + [-axisSize/2 axisSize/2]]);
    axis (imax,'equal'); 
    axis (imax,[loc(1,j) + [-axisSize/2 axisSize/2], loc(2,j) + [-axisSize/2 axisSize/2]]);
    hold (imax,'off');
    set(imax, 'XTick', [], 'YTick', [], 'Color', 'k');
    if (~isempty(track.run))
        xl = get(imax, 'XLim');
        xl = xl(1)*0.99 + xl(2)*0.01;
        yl = get(imax, 'YLim');
        yl = yl(1)*0.01 + yl(2)*0.99;
        
        t = {['\color{orange}time = ' num2str(eti(sind(j)), '%.1f')]};
        I = find([track.run.startInd] <= sind(j) & [track.run.endInd] >= sind(j), 1, 'first');
        %col = {[0.6 0.6 .9], [1 1 1], [0.7 1 0.7],[1 0.7 0.7]};
    
         if (~isempty(I))
            t = [t ['\color[rgb]{0.6 0.6 0.9}run, \tau = ' num2str(track.run(I).runTime,'%.1f') ' s, L = ' num2str(track.run(I).pathLength, '%.1f') ' cm']];
        end
       
        I = find([track.reorientation.startInd] <= sind(j) & [track.reorientation.endInd] >= sind(j), 1, 'first');
        if (~isempty(I))
            t = [t ['\color{white}turn, \Delta\theta = ' num2str(rad2deg(diff(unwrap([track.reorientation(I).prevDir;track.reorientation(I).nextDir]))), '%.1f')]];
        end
        
        I = find([track.headSwing.startInd] <= sind(j) & [track.headSwing.endInd] >= sind(j), 1, 'first');
        if (~isempty(I))
            if (track.headSwing(I).accepted)
                t = [t '\color[rgb]{0.7 1 0.7}accepted head sweep'];
            else
                t = [t '\color[rgb]{1 0.7 0.7}rejected head sweep'];
            end
        end
    end
    if (~isempty(t))
        text(double(xl), double(yl), t, 'Interpreter', 'Tex', 'Parent', imax, 'FontName', 'Arial', 'FontSize', fontsize, 'HorizontalAlignment', 'Left', 'VerticalAlignment', 'Top', 'BackgroundColor', 'k');
    end
    
    handles = updateCenter(handles, track, sind(j), sstart, send, dataaxes);
    
    
    if (~isempty(aviobj) || ~isempty(vidObj))
        if (isempty(avirect))
            F = getframe(gcf);
        else
            F = getframe(gcf, avirect);
        end        
        if (~isempty(aviobj))
            aviobj = addframe(aviobj, F);
        end
        if (~isempty(vidObj))
            writeVideo(vidObj,F);
        end
    else
        timeleft = delayTime - toc(ts1);
        if (timeleft > 0)
            pause(timeleft);
        else
            pause(0.001);
        end
    end
        
end
%{
drawtime
tracktime
otherplotstime
%}
end

function datafields(track, inds, dataaxes, fontsize)
    thetafield = 'sspineTheta';
%     if (false && strcmpi (track.so.method,'new'))
%         fields = {thetafield,track.so.speed_field,'spheadperp'};
%         ftitles = {'Bondy Bend Angle', track.so.speed_field, 'vhead perp'};
%         mult = [rad2deg(1), 1, 1];
%     else
%         fields = {thetafield,track.so.speed_field,'vel_dp'};
%         ftitles = {'Bondy Bend Angle', track.so.speed_field, 'Velocity Dot Product'};
%         mult = [rad2deg(1), 1, 1];
%     end

    fields = {thetafield, 'vel_dp', track.so.speed_field };
    ftitles = {'$\theta body $', '$v{\cdot}h$','cm/min'};
    mult = [180/pi 1 60];
    
    
    
    start{1} = [track.run.startInd];
    stop{1} = [track.run.endInd];
    start{2} = [track.reorientation([track.reorientation.numHS] >= 0).startInd];
    stop{2} = [track.reorientation([track.reorientation.numHS] >= 0).endInd];
    start{3} = [track.headSwing([track.headSwing.accepted]).startInd];
    stop{3} = [track.headSwing([track.headSwing.accepted]).endInd];
    start{4} = [track.headSwing(~[track.headSwing.accepted]).startInd];
    stop{4} = [track.headSwing(~[track.headSwing.accepted]).endInd];

    %fn = {'run', 'turn', 'ahs', 'rhs'};
    col = {[0.6 0.6 .9], [1 1 1], [0.7 1 0.7],[1 0.7 0.7]};
    
    for j = 1:length(start)
        start{j} = start{j}(start{j} > inds(1));
        stop{j} = stop{j}(stop{j} > start{j}(1) & stop{j} < inds(end));
        start{j} = start{j}(start{j} < stop{j}(end));
    end
    
    eti = track.getDerivedQuantity('eti');
    for k = 1:3
        for m = 1:length(start)
            yl = [min(mult(k)*track.dq.(fields{k})(inds)) max(mult(k)*track.dq.(fields{k})(inds))];
            dy = diff(yl)/25;
            
            si = eti(start{m});
            ei = eti(stop{m});
            c = col{m};
            for j = 1:length(si)
                patch([si(j) ei(j) ei(j) si(j) si(j)], [yl(1) yl(1) yl(2) yl(2) yl(1)], c, 'EdgeColor', 'none', 'Parent', dataaxes(k)); hold(dataaxes(k), 'on');
            end
        end

        
        plot (dataaxes(k), track.dq.eti(inds), mult(k)*track.dq.(fields{k})(inds), 'k-', 'LineWidth', 2); hold(dataaxes(k), 'on');
        
        ylabel(dataaxes(k), ftitles{k}, 'Interpreter', 'Latex', 'FontSize', fontsize, 'FontName', 'Arial');
    end
    
    spfields = {{'theta_cut', 'headswing_start', 'headswing_stop'}, {'aligned_dp'},{'stop_speed_cut', 'start_speed_cut'}};
    spcolors = {{'m-','g-','r-'}, {'m-'},{'r-', 'g-'}};
    spmirror = [1 0 0];
    
    for k = 1:3
        
        x = track.dq.eti(inds);
        for j = 1:length(spfields{k})
            f = spfields{k}{j};
            c = spcolors{k}{j};
            y = repmat (mult(k)*track.so.(f), size(x));
            plot (dataaxes(k), x,y,c);
            if (spmirror(k));
                plot (dataaxes(k), x,-y,c);
            end
            yl = [min(mult(k)*track.dq.(fields{k})(inds)) max(mult(k)*track.dq.(fields{k})(inds))];
            dy = diff(yl)/25;
            ylim(dataaxes(k), yl + dy*[-1 1]);
            set(dataaxes(k), 'box', 'on');
        end
    end
    p = get(dataaxes(1), 'position');
    py = -0.95*p(2)/p(4);
    hx = xlabel(dataaxes(1), 'time (s)');
    set(hx, 'Units', 'normalized');
    p = get(hx, 'position');
    p(2) = py;
    set (hx, 'position', p, 'VerticalAlignment', 'bottom');
end

function handles = updateCenter(handles, track, cind, start, stop, dataaxes)
    thetafield = 'sspineTheta';
    if ~isempty(handles)
        for j = 1:length(handles)
            delete(handles(j));
        end
    end
    fields = {thetafield, 'vel_dp', track.so.speed_field };
   % ftitles = {'$\theta_{body}$', '$\hat{v}\cdot\hat{h}$','speed (cm/min)' };
    mult = [180/pi 1 60];
   
    for k = 1:3
        %subplot(2,2,k+1); hold on
        hold(dataaxes(k), 'on');
        handles(k) = plot (dataaxes(k), track.dq.eti(cind), mult(k)*track.dq.(fields{k})(cind), 'c.','MarkerSize',25);
        xlim(dataaxes(k), [min(track.dq.eti(start)) max(track.dq.eti(stop))]);        
    end
end   

function [imax, dataaxes] = setupFigure(fignum, resizefig)
    existsAndDefault('fignum', gcf);
    clf(fignum);
    if (resizefig)
        p = get(fignum, 'position');
        pold = p;
        p(3) = p(4) * 8/3;
        p(1) = max(1/2*(pold(3)-p(3)) + p(1), 0);
        set(fignum, 'position', p);
    end
    h = 0.8;
    w = h * 3/8;
    imax = axes('position', [(.5-w)/2 (1-h)/2 w h], 'Parent', fignum);
    
    h = 0.9;
    w = 0.45;
    for j = 1:3
        dataaxes(j) = axes('position', [.5+(.5-w)/2 0.9*(1-h)+(j-1)*h/3 w h/3.1], 'Parent', fignum); %#ok<AGROW>
    end
    for j = 2:3
        set(dataaxes(j), 'XTickLabel', []);
    end
end



