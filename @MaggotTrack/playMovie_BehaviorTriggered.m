function playMovie_BehaviorTriggered(track, fieldName,fieldDescription, displacementAxis, triggeredSum, triggeredInd, figTitle, varargin)
%@MaggotTrack
%playMovie(track, varargin)
%enter options as pairs, caps matter
%options, with defaults
%
%frameRate = [];
%ptbuffer = 1000;
%delayTime = 0.05;
%axisSize (3 * length of larva)
%inds = 1:length(track.pt);
%startLoc = []; > if startLoc & stopLoc are not empty, we run the movie
%between these two points
%stopLoc = []; >
%startTime = []; if startTime and stopTime are not empty, re run the movie
%between these two times
%stopTime = [];
%vidObj = [];
%aviobj = [];
%avirect = [];
%resizefig = true;
%figcolor = 'w';
%fontcolor = 'k';
%underlayTemporalField = ''; boolean true if on
%underlayTemporalOnColor = [0 0.8 0.8];
%underlayTemporalOffColor = [0 0 0];
%underlayTemporalOnMessage = 'Light on';
%underlayTemporalOffMessage = 'Light off';
%rect = [];
%mark_size = 5

%pass 'fid', [] to not load images from disk
%ptskip, [] -- if not empty we sample pts at this interval 
if (length(track) > 1)
    for j = 1:length(track)
        track(j).playMovie(varargin{:});
    end
    return;
end
frameRate = []; %override frameRate to set a target frame rate
stretchReos = false; %if skipping frames due to frame rate, then don't skip reorientation frames

ptbuffer = ceil(40/track.dr.interpTime);
delayTime = max(track.dr.interpTime/4-0.1, 0.01);
%set(0,'DefaultTextInterpreter', 'Latex');
if (isempty(track.expt) || isempty(track.expt.camcalinfo))
    mf = 1;
else
    mf = track.expt.camcalinfo.realUnitsPerPixel;
end
%axisSize = mf*max(size(track.pt(1).imData));
sl = track.getDerivedQuantity('spineLength');
axisSize = 3*median(sl(:,track.getDerivedQuantity('ihtValid')));
if (axisSize <= 0)
    try 
        track.expt.openDataFile();
        pttemp = track.expt.reloadPoint(track.pt(1));
        axisSize = mf*max(size(pttemp.imData));
    catch
    end
end


if (axisSize <= 0 || ~isfinite(axisSize))
    %axisSize = 50;
    axisSize = 8 * mean(track.getDerivedQuantity('speed'));
     
end
ptskip = [];
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
scalebar = false;
mark_size = 5;

if (~isempty(track.expt))
    track.expt.openDataFile;
    fid = track.expt.fid;
else
    fid = [];
end
if (isa(track.pt(1), 'MWTTrackPoint'))
    imOptions = {'pretty', true, 'drawHeadArrow', true, 'spineColor', 'y.','contourColor', 'm-', 'contourWidth', 1, 'scale', 1, 'drawContour', true, 'drawSpine', true};%'drawContour', true, 'scale', 1,'contourColor', 'r-', 'LineWidth', 1, 'mhWidth', 1};
else
    imOptions = {'pretty', true, 'drawHeadArrow', false, 'spineColor', 'y.','contourColor', 'c.', 'contourWidth', 1, 'scale', 1, 'drawContour', false, 'drawSpine', true};%'drawContour', true, 'scale', 1,'contourColor', 'r-', 'LineWidth', 1, 'mhWidth', 1};
end

underlayImData = [];
overlayImData = [];
underlayTemporalField = '';
underlayTemporalOnColor = [0 1 1];
underlayTemporalOffColor = [0 0 0];
underlayTemporalOnMessage = 'Light on';
underlayTemporalOffMessage = 'Light off';
rect = [];
fixedrect = false;
varargin = assignApplicable(varargin);
if (~isempty(iinds))
    inds = track.getDerivedQuantity('mapinterpedtopts', false,iinds);
end
if (~isempty(ptskip))
    inds = inds(1):ptskip:inds(end);
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

    
%pt = [track.pt];
sstart = sind(1) - ptbuffer;
send = sind(end) + ptbuffer;
if (sstart < 1)
    sstart = 1;
end
if (send > length(sloc))
    send = length(sloc);
end

%dounderlay = ~isempty(underlayImData);
dounderlay = false;
[imax, dataaxes, underlayax] = setupFigure(gcf, resizefig,dounderlay);

u = get(imax, 'units');
set(imax, 'units', 'pixels');
imh = get(imax, 'position');
imh = imh(4);
fontsize = .04 * imh;
set(imax, 'units', u);
%set(gcf, 'DefaultAxesColor', 'w', 'DefaultTextColor', 'k');
set(gcf, 'Color', figcolor);
%set(imax, 'Color', 'k');
set(dataaxes, 'XColor', fontcolor, 'YColor', fontcolor, 'FontSize', fontsize-2)
datafields(track, sstart:send, dataaxes, fontsize, fieldName, fieldDescription, displacementAxis, triggeredSum);

% if (dounderlay)
%    % set(pcolor (underlayax, underlayImData.x, underlayImData.y, underlayImData.im), 'FaceAlpha', 0.4, 'EdgeColor', 'none');
%     imagesc(underlayImData.x, underlayImData.y, double(underlayImData.im), 'Parent', underlayax);
%     imOptions = [imOptions {'patchOptions', {'AlphaData', 0.4}}];
%     axis(underlayax, 'equal'); axis(underlayax, 'xy'); axis(underlayax, 'off');
% else
%     overlayax = [];
% end

handles = [];
if (isempty(track.expt))
    ccinfo = [];
else
    ccinfo = track.expt.camcalinfo;
end
eti = track.getDerivedQuantity('eti');
temporalunderlay = false;
if (~isempty(underlayTemporalField))
    temporalFieldData = logical(setNonFiniteToZero(track.getDerivedQuantity(underlayTemporalField)));
    temporalunderlay = true;
end

if (fixedrect && isempty(rect))
    plot (imax, loc(1,inds) - axisSize/2, loc(2,inds) - axisSize/2, loc(1,inds) + axisSize/2, loc(2,inds) + axisSize/2);
    axis(imax, 'equal');
    rect = [get(imax, 'xlim') get(imax, 'ylim')];
    ll = min(loc(:,inds),[],2);% - axisSize/2;
    ur = max(loc(:,inds),[],2);% + axisSize/2;
   
    if (ll(1) < rect(1) || ur(1) > rect(2))
        rect(1:2) = rect(1:2) + mean([ll(1) ur(1)]) - mean(rect(1:2));
    end
    if (ll(2) < rect(3) || ur(2) > rect(4))
        rect(3:4) = rect(3:4) + mean([ll(2) ur(2)]) - mean(rect(3:4));
    end
    cla(imax);
    
end

if (~isempty(underlayImData))
    if (isempty(rect))
        ll = min(loc(:,inds),[],2);% - axisSize/4;
        ur = max(loc(:,inds),[],2);% + axisSize/4;
    else
        ll = [rect(1) rect(3)];
        ur = [rect(2) rect(4)];
    end
    xv = find(underlayImData.x >= ll(1) & underlayImData.x <= ur(1));
    yv = find(underlayImData.y >= ll(2) & underlayImData.y <= ur(2));
   % [xx,yy] = meshgrid(xv,yv);
    underlayImData.x = underlayImData.x(xv);
    underlayImData.y = underlayImData.y(yv);
    underlayImData.im = underlayImData.im(yv,xv,:);
    imOptions = [imOptions, {'underlayImData', underlayImData}];
    size(underlayImData.im)
end

ts0 = tic();
nframes = 0;
skipCount = 1;
%if (~isempty(frameRate) && stretchReos)
 %   disp ('stretching reorientations out to slow speed; runs will animate faster');
%end
axes(imax); axis xy;
triggeredInd = track.getDerivedQuantity('mapinterpedtopts', false, triggeredInd);
triggeredInd = sort(triggeredInd, 'ascend');
triggeredInd = triggeredInd(triggeredInd > inds(1));
for j = inds
    nframes = nframes + 1;
    I = find([track.reorientation.startInd] <= sind(j) & [track.reorientation.endInd] >= sind(j), 1, 'first');
    if ((mod(nframes, skipCount) >= 1) && ~(stretchReos && ~isempty(I)))
        continue;
    end
    ts1 = tic();
    hold (imax,'off'); 
    cla(imax); 
   % set(imax, 'Color', 'k');
   if (isempty(rect))
       imrect = [loc(1,j) + [-axisSize/2 axisSize/2], loc(2,j) + [-axisSize/2 axisSize/2]];
   else
       imrect = rect;
   end
   if (temporalunderlay)
       if (~isempty(ccinfo))
           uid.x = imrect(1):ccinfo.realUnitsPerPixel:imrect(2);
           uid.y = imrect(3):ccinfo.realUnitsPerPixel:imrect(4);           
       else
           uid.x = round(imrect(1)):round(imrect(2));
           uid.y = round(imrect(3)):round(imrect(4));
       end
       uid.im = ones([length(uid.y) length(uid.x) 3]);
       for k = 1:3
           if (temporalFieldData(sind(j)))
               uid.im(:,:,k) = underlayTemporalOnColor(k);
           else
               uid.im(:,:,k) = underlayTemporalOffColor(k);
           end
       end
       pt(j).drawTrackImage(ccinfo,'fid', fid, 'Axes', imax, imOptions{:}, 'underlayImData', uid, 'underlayScale', 0.6);
   else
       pt(j).drawTrackImage(ccinfo,'fid', fid, 'Axes', imax, imOptions{:}); 
       set(imax, 'YDir', 'normal');
   end
    hold (imax, 'on');
 
    sstart = sind(j) - ptbuffer;
    send = sind(j) + ptbuffer;
    if (sstart < 1)
        sstart = 1;
    end
    if (send > length(sloc))
        send = length(sloc);
    end
    
    
        
    
    plot (imax,sloc(1,sstart:send), sloc(2,sstart:send), 'w.-', 'MarkerSize', mark_size);
%    plot (imax, sloc(1, sstart:(1/track.dr.interpTime):send), sloc(2, sstart:(1/track.dr.interpTime):send), 'w.');
    plot (imax,sloc(1,sind(j)), sloc(2,sind(j)), 'wo', 'MarkerSize', mark_size);
    
    axis (imax,imrect);
    axis (imax,'equal'); 
    axis (imax,imrect);
    
    set(imax, 'XTick', [], 'YTick', []);%, 'Color', 'k');
    if (scalebar && ~isempty(ccinfo))
        sbsize = double(ceil(axisSize * 2.5));
        xl = get(imax, 'XLim');
        yl = get(imax, 'YLim');
        sbx = xl(2) + [-1.1 -0.1]*sbsize/10;
        sby = yl(1)*0.95 + yl(2)*0.05;
        plot (imax, sbx, [sby sby], 'w-', 'LineWidth', 4);
        text (mean(sbx), yl(1)*0.94 + yl(2)*0.06, [num2str(sbsize) ' mm'], 'Color', 'w',...
            'Interpreter', 'Tex', 'Parent', imax, 'FontName', 'Arial', 'FontSize', fontsize, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom',...
            'BackgroundColor', 'k');
    end
    hold (imax,'off');
    xl = get(imax, 'XLim');
    xl = xl(1)*0.99 + xl(2)*0.01;
    yl = get(imax, 'YLim');
    yl = yl(1)*0.01 + yl(2)*0.99;
    t = {['\color{orange}time = ' num2str(eti(sind(j)), '%.1f')]};
    if (temporalunderlay)
        if (temporalFieldData(sind(j)))
            t = [t '\color{white}' underlayTemporalOnMessage];
        else
            t = [t '\color{white}' underlayTemporalOffMessage];
        end
    end
    if (~isempty(track.run))
        
        
        I = find([track.run.startInd] <= sind(j) & [track.run.endInd] >= sind(j), 1, 'first');
        %col = {[0.6 0.6 .9], [1 1 1], [0.7 1 0.7],[1 0.7 0.7]};
    
         if (~isempty(I))
            t = [t ['\color[rgb]{0.6 0.6 0.9}run, \tau = ' num2str(track.run(I).runTime,'%.1f') ' s, L = ' num2str(track.run(I).pathLength, '%.1f') ' cm']];
        end
       
        I = find([track.reorientation.startInd] <= sind(j) & [track.reorientation.endInd] >= sind(j), 1, 'first');
        if (~isempty(I))
            if (track.reorientation(I).numHS == 0)
                t = [t ['\color{white}pause, \Delta\theta = ' num2str(rad2deg(diff(unwrap([track.reorientation(I).prevDir;track.reorientation(I).nextDir]))), '%.1f')]];
            else
                t = [t ['\color{white}turn, \Delta\theta = ' num2str(rad2deg(diff(unwrap([track.reorientation(I).prevDir;track.reorientation(I).nextDir]))), '%.1f')]];
            end
        end
        
        I = find([track.headSwing.startInd] <= sind(j) & [track.headSwing.endInd] >= sind(j), 1, 'first');
        if (~isempty(I))
            if (~track.headSwing(I).valid)
                inv = ' *invalid';
            else
                inv = '';
            end
            if (track.headSwing(I).accepted)
                t = [t '\color[rgb]{0.7 1 0.7}accepted head sweep' inv];
            else
                t = [t '\color[rgb]{1 0.7 0.7}rejected head sweep' inv];
            end
        end
    end
    if (~isempty(t))
        text(xl, yl, t, 'Interpreter', 'Tex', 'Parent', imax, 'FontName', 'Arial', 'FontSize', fontsize, 'HorizontalAlignment', 'Left', 'VerticalAlignment', 'Top', 'BackgroundColor', 'k');
    end
    
    handles = updateCenter(handles, track, sind(j), find(track.dq.eti < track.dq.eti(sind(j)) + displacementAxis(1), 1, 'last'), find(track.dq.eti > track.dq.eti(sind(j)) + displacementAxis(end), 1, 'first'), dataaxes);
    if (~isempty(underlayax))
        set(underlayax, 'Position', get(imax, 'position'), 'XLim', get(imax, 'XLim'), 'YLim', get(imax, 'YLim')); 
    end
    
    if (dounderlay)
        set(imax, 'Color', 'none');
    else
        set(imax, 'Color', 'k');
    end
    if (~isempty(figTitle))
        title(imax, figTitle,'FontName', 'Arial', 'FontSize', fontsize);
    end
    triggeredInd = triggeredAnimation(track, triggeredInd, j, dataaxes, fieldName, displacementAxis, triggeredSum, aviobj, vidObj);
    
    
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
        if (~isempty(frameRate))
            vidTime = toc(ts0);
            if (vidTime > 2) %update frame rate every 2 seconds
                fr = nframes/vidTime;
                nframes = 0;
                ts0 = tic();
                if (fr < frameRate)
                    delayTime = 0.001;
                    skipCount = skipCount * min(2,frameRate/fr);
                else
                    skipCount = max(1,skipCount * frameRate/fr);
                    if (skipCount == 1)
                        delayTime = 1/frameRate;
                    end
                end
            end
        end
    end
    
    
        
end
%{
drawtime
tracktime
otherplotstime
%}
end

function datafields(track, inds, dataaxes, fontsize, fieldName, fieldDescription, displacementAxis, triggeredSum)
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

    fields = {thetafield, track.so.speed_field, fieldName};
    ftitles = {'$\theta body $',{'speed', 'cm/min'}, fieldDescription};
    mult = [180/pi 60 1];
    
    
    if (isempty(track.run) || isempty(track.reorientation) || isempty(track.headSwing))
        start = {[],[],[],[]};
        stop = {[],[],[],[]};
    else

        start{1} = [track.run.startInd];
        stop{1} = [track.run.endInd];
        start{2} = [track.reorientation([track.reorientation.numHS] >= 0).startInd];
        stop{2} = [track.reorientation([track.reorientation.numHS] >= 0).endInd];
        start{3} = [track.headSwing([track.headSwing.accepted]).startInd];
        stop{3} = [track.headSwing([track.headSwing.accepted]).endInd];
        start{4} = [track.headSwing(~[track.headSwing.accepted]).startInd];
        stop{4} = [track.headSwing(~[track.headSwing.accepted]).endInd];
    end

    %fn = {'run', 'turn', 'ahs', 'rhs'};
    col = {[0.6 0.6 .9], [1 1 1], [0.7 1 0.7],[1 0.7 0.7]};
    
    for j = 1:length(start)
        if (isempty(start{j}) || isempty(stop{j}))
            continue;
        end
        start{j} = start{j}(start{j} > inds(1));
        stop{j} = stop{j}(stop{j} > start{j}(1) & stop{j} < inds(end));
        start{j} = start{j}(start{j} < stop{j}(end));
    end
    
    eti = track.getDerivedQuantity('eti');
    
    plot (dataaxes(1), displacementAxis, triggeredSum, 'k-', 'LineWidth', 2);
    set (dataaxes(1), 'YTick', []);
    ylabel(dataaxes(1), ['$\langle$' fieldDescription '$\rangle$'], 'Interpreter', 'Latex', 'FontSize', fontsize, 'FontName', 'Arial');
    yl = get(dataaxes(1), 'YLim'); ylim(dataaxes(1), [-1 1] * max(abs(yl)));
    for k = 1:2
        for m = 1:length(start)
            yl = [min(mult(k)*track.dq.(fields{k})(inds)) max(mult(k)*track.dq.(fields{k})(inds))];
           
            
            si = eti(start{m});
            ei = eti(stop{m});
            c = col{m};
            for j = 1:length(si)
                patch([si(j) ei(j) ei(j) si(j) si(j)], [yl(1) yl(1) yl(2) yl(2) yl(1)], c, 'EdgeColor', 'none', 'Parent', dataaxes(k+1)); hold(dataaxes(k+1), 'on');
            end
        end

        
        plot (dataaxes(k+1), track.dq.eti(inds), mult(k)*track.dq.(fields{k})(inds), 'k-', 'LineWidth', 2); hold(dataaxes(k), 'on');
        
        ylabel(dataaxes(k+1), ftitles{k}, 'Interpreter', 'Latex', 'FontSize', fontsize, 'FontName', 'Arial');
    end
    
    k = 3;
     plot (dataaxes(k+1), track.dq.eti(inds), mult(k)*track.dq.(fields{k})(inds), 'k-', 'LineWidth', 2); hold(dataaxes(k), 'on');
        
        ylabel(dataaxes(k+1), ftitles{k}, 'Interpreter', 'Latex', 'FontSize', fontsize, 'FontName', 'Arial');
        
    set(dataaxes(4), 'YTick', [], 'XTickLabel', []);
    
    spfields = {{'theta_cut', 'headswing_start', 'headswing_stop'},{'stop_speed_cut', 'start_speed_cut'}, {}};
    spcolors = {{'m-','g-','r-'},{'r-', 'g-'}, {}};
    spmirror = [1 0 0];
    
    for k = 1:3
        
        x = track.dq.eti(inds);
        for j = 1:length(spfields{k})
            f = spfields{k}{j};
            c = spcolors{k}{j};
            y = repmat (mult(k)*track.so.(f), size(x));
            plot (dataaxes(k+1), x,y,c);
            if (spmirror(k));
                plot (dataaxes(k+1), x,-y,c);
            end
            yl = [min(mult(k)*track.dq.(fields{k})(inds)) max(mult(k)*track.dq.(fields{k})(inds))];
            dy = diff(yl)/25;
            ylim(dataaxes(k+1), yl + dy*[-1 1]);
            set(dataaxes(k+1), 'box', 'on');
        end
        if (k == 3)
            ylim(dataaxes(k+1), [-1 1]*percentile(abs(track.dq.(fields{k})), .99));
        end
    end
    p = get(dataaxes(2), 'position');
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
    
    fields = {thetafield, track.so.speed_field};
    mult = [180/pi 60];
   
    for k = 1:3
        %subplot(2,2,k+1); hold on
        hold(dataaxes(k+1), 'on');
        if (k == 3)
            handles(k) = plot (dataaxes(k+1), track.dq.eti(cind)*[1 1], get(dataaxes(k+1), 'YLim'), 'c--','LineWidth',4);
        else
            handles(k) = plot (dataaxes(k+1), track.dq.eti(cind), mult(k)*track.dq.(fields{k})(cind), 'c.','MarkerSize',25);
        end
        if (isempty(start))
            start = 1;
        end
        if (isempty(stop))
            stop = length(track.dq.eti);
        end
        xlim(dataaxes(k+1), [min(track.dq.eti(start)) max(track.dq.eti(stop))]);        
    end
    
end   

function [imax, dataaxes, underlayax] = setupFigure(fignum, resizefig, makeunderlay)
    existsAndDefault('fignum', gcf);
    existsAndDefault('makeunderlay', false);
    existsAndDefault('resizefig', true);
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
    if (makeunderlay)
        underlayax = axes('position', [(.5-w)/2 (1-h)/2 w h], 'Parent', fignum, 'XTick', [], 'YTick', [], 'Color', 'k');
        imax = axes('position', [(.5-w)/2 (1-h)/2 w h], 'Parent', fignum, 'Color', 'none', 'XTick', [], 'YTick', []);
    else
        imax = axes('position', [(.5-w)/2 (1-h)/2 w h], 'Parent', fignum, 'XTick', [], 'YTick', [], 'Color', 'k');
        underlayax = [];
    end
    h = 0.9;
    w = 0.45;
    for j = 1:4
        dataaxes(j) = axes('position', [.5+(.5-w)/2 0.9*(1-h)+(j-1)*h/4 w h/4.1], 'Parent', fignum); %#ok<LAXES,AGROW>
    end
    for j = 2:4
        set(dataaxes(j), 'XTickLabel', []);
    end
end

function triggeredInd = triggeredAnimation(track, triggeredInd, cind, dataaxes, fieldName, displacementAxis, triggeredSum, aviobj, vidObj)
    persistent triggered bc1 c0 p0 xl0 xd0;
    maxcount = 30;
   
    if (isempty(triggeredInd) || cind < min(triggeredInd))
        triggered = false;
        return;
    end
    
    count = cind - min(triggeredInd);
    if (count > maxcount)
        c = get(dataaxes(end), 'Children'); 
        set(dataaxes(end), 'Position', p0, 'Color', bc1);
        set(c(end), 'Color', c0);
        triggeredInd = triggeredInd(triggeredInd > cind);
        triggered = false;
        return;
    end
    
    
    
    transitioncount = 20;
    c = get(dataaxes(end), 'Children'); 
    if (~triggered)
         bc1 = get(dataaxes(end), 'Color');   
         c0 = get(c(end), 'Color');
         xd0 = get(c(1), 'XData');
         p0 = get(dataaxes(end), 'position');
         xl0 = get(dataaxes(end), 'XLim');
         triggered = true;
    end
    set(c(1), 'XData', xd0);
    set(dataaxes(end), 'Color', 'none','XTickLabel', []);
    set(c(end), 'Color', 'r');
    pf = get(dataaxes(1), 'position');
    cc = min(count, transitioncount);
    ind0 = find(track.getDerivedQuantity('eti') >= track.pt(cind-count).et,1,'first');%track.getDerivedQuantity('mapptstointerped', false, cind-count);
    
    set(dataaxes(end), 'Position', (cc*pf + (transitioncount - cc)*p0)/transitioncount, 'XLim', track.dq.eti(ind0) + [min(displacementAxis) max(displacementAxis)]);%, 'YLim', (j*ylf + (nsteps - j)*yl0)/nsteps);
    if (transitioncount == count) 
            set(get(dataaxes(1), 'Children'), 'YData', triggeredSum + interp1(track.dq.eti - track.dq.eti(ind0), track.dq.(fieldName), displacementAxis, 'linear', 0));
    end
    %{
    else
            set(dataaxes(end), 'Position', p0, 'Color', bc1, 'YLim', yl0);
            set(c(end), 'Color', c0);
            set(get(dataaxes(1), 'Children'), 'YData', triggeredSum + interp1(track.dq.eti - track.dq.eti(cind), track.dq.(fieldName), displacementAxis, 'linear', 0));
        end
%}
        %{
        for kk = 1:nreps
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
                pause(1/30);
            end
        end
        %}
  %  end

end