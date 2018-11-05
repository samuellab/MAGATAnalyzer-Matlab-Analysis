function  toMWTFiles(expt, prefix, makeMultipleBlobsFile, varargin) %outputs files to be used in mwt choreography analysis
%function  toMWTFiles(expt, prefix, makeMultipleBlobsFile, varargin) %outputs files to be used in mwt choreography analysis
%creates MWT files (.summary and .blob or .blobs that can be loaded into
%                   MWT choreography)
%
%expt < Experiment -- must have MaggotTracks
%prefix -- the prefix to be used; 
%          default:  expt.fname = dirToExpt/fname.bin
%                then prefix = dirToExpt/MWT/fname
%makeMultipleBlobsFile -- if true, we package all the tracks into a single
%                                  .blobs file
%                         if false(default), each track is in its own .blob
%                                   file
%
%{
    MWT SUMMARY FILE : our approximation
    1. Image number (counts up from 1) : Same
    2. Image time (in seconds from the beginning of the experiment) : Same
    3. Number of objects tracked : Same
    4. Number of objects that have lasted long enough to be averaged into summary data : Same as 3
    5. Average duration for an object to have been tracked : min(time, average track duration)
    6. Average speed in pixels/second : mean(speed(t))
    7. Average angular speed in radians/second : mean abs(dtheta/dt (t))
    8. Average length of object in pixels : mean (spineLength(t))
    9. Average relative length of object, where relative length = current length / mean length :  8/mean(8)
    10. Average width of object in pixels : area / spinelength(t)
    11. Average relative width of object : 10 / mean(10)
    12. Average aspect ratio of object (width / length) : 10/8
    13. Average relative aspect ratio (current W/L) / (mean W/L) : 12 / mean(12)
    14. Average “end wiggle”: angle in radians between the last 20% of the body and the rest of the
    body (using whichever end shows a greater angle). : mean(sspineTheta(t))
    15. Average number of pixels filled as part of object detection : mean(area(t))

    %% followed by pairs of numbers that specify how objects were found and lost. The first
    number indicates the object number before this image was processed; the second, what it
    turned into afterwards. 0 means that no object was found. Thus, 0 24 means a new object
    was found and was given number 24; 24 0 means that the existing object 24 was lost; 24 40
    25 40 means that 24 and 25 collided and the resulting object was labeled 40; and 40 67 68
    means that compound object 40 fell apart into two objects which were labeled 67 and 68.

    %%% followed by pairs of numbers if and only if objects were grouped into blobs files.
        The numbers are in the form: N X.Y, where N is the object number, X is the file number,
        and Y is the byte offset within the file.
%}

giveUpdates = true;
varargin = assignApplicable(varargin);

existsAndDefault('makeMultipleBlobsFile', false);
[d,f] = fileparts(expt.fname);
existsAndDefault ('prefix', fullfile(d, 'MWT', f));

if (~exist(fileparts(prefix), 'dir'))
    mkdir(fileparts(prefix));
end

startFrame = [expt.track.startFrame] + 1;
endFrame = [expt.track.endFrame] + 1;
trackNum = [expt.track.trackNum];
npts = [expt.track.npts];

filenum = zeros(size(endFrame));
locInFile = filenum;
ts1 = tic; lastt = toc(ts1);
if(makeMultipleBlobsFile)
    fid = fopen([prefix '_00000k.blobs'], 'wt');
    [~,I] = sort(endFrame);
    for j = 1:length(I)
        if (mod(j,1000) == 0)
            fclose(fid);
            fid = fopen([prefix sprintf('_%05dk.blobs',floor(j/1000))], 'wt');
        end
        tn = I(j);        
        filenum(tn) = floor(j/1000);
        locInFile(tn) = ftell(fid);
        fprintf(fid, '%s', expt.track(j).toMWTTrackString);
        if (giveUpdates && toc(ts1) - lastt > 60 )
            lastt = toc(ts1);
            fprintf('%.0f sec elapsed; %d / %d tracks written; %.0f %% of pts written; estimated time remaining = %.1f min\n', ...
                lastt, j, length(I), 100*sum(npts(I(1:j)))/sum(npts), lastt * (sum(npts) - sum(npts(I(1:j))))/sum(npts(I(1:j)))/60);
        end
    end
    fclose(fid);
else
    for j = 1:length(expt.track)
        expt.track(j).toMWTTrackFile(prefix);
        if (giveUpdates && toc(ts1) - lastt > 60 )
            lastt = toc(ts1);
            fprintf('%.0f sec elapsed; %d / %d tracks written; %.0f %% of pts written; estimated time remaining = %.1f min\n', ...
                lastt, j, length(expt.track), 100*sum(npts(1:j))/sum(npts), lastt * (sum(npts) - sum(npts(1:j)))/sum(npts(1:j))/60);
        end
    end
end
        
%make summary file        
if (giveUpdates)
    disp([num2str(toc(ts1)), ' making summary file']);
end
if (~makeMultipleBlobsFile)
    fileNum = [];
    locInFile = [];
end
makeMWTSummaryFile(expt, prefix, fileNum, locInFile, 'giveUpdates', giveUpdates, varargin{:});
% et = expt.elapsedTime;
% frameNum = 1:length(et);
% numTracked = zeros(size(frameNum));
% 
% for j = 1:length(startFrame)
%     numTracked(startFrame(j):endFrame(j)) = numTracked(startFrame(j):endFrame(j)) + 1;
% end
% 
% 
% trackDuration = et(endFrame) - et(startFrame);
% dt = median(diff(et));
% binEdges = [et(1)-dt/2, 0.5*(et(1:end-1)+et(2:end)) et(end)+dt/2];
% 
% %convert lengths to pixels if not already there
% if ~isempty(expt.camcalinfo)
%     pm = expt.camcalinfo.pixelsPerRealUnit;
% else
%     pm = 1;
% end
% 
% eti = expt.gatherField('eti');
% sp = pm*expt.gatherField('speed');
% sl = pm*expt.gatherField('spineLength');
% dtheta = expt.gatherField('deltatheta');
% ia = expt.gatherField('iarea');
% width = ia./sl;
% st = expt.gatherField('sspineTheta');
% scr = expt.gatherField('scovRatio');
% 
% [~,spvstime] = meanyvsx(eti(isfinite(sp)), sp(isfinite(sp)), binEdges);
% [~,slvstime] = meanyvsx(eti(isfinite(sl)), sl(isfinite(sl)), binEdges);
% [~,dthetavstime] = meanyvsx(eti(isfinite(dtheta)), abs(dtheta(isfinite(dtheta))), binEdges);
% [~,areavstime] = meanyvsx(eti(isfinite(ia)), ia(isfinite(ia)), binEdges);
% [~,widthvstime] = meanyvsx(eti(isfinite(width)), width(isfinite(width)),binEdges);
% [~,arvstime] = meanyvsx(eti(isfinite(scr)), 1./scr(isfinite(scr)), binEdges);
% [~,stvstime] = meanyvsx(eti(isfinite(st)), abs(st(isfinite(st))), binEdges);
% 
% %{
% MWT SUMMARY FILE : our approximation
%     1. Image number (counts up from 1) : Same
%     2. Image time (in seconds from the beginning of the experiment) : Same
%     3. Number of objects tracked : Same
%     4. Number of objects that have lasted long enough to be averaged into summary data : Same as 3
%     5. Average duration for an object to have been tracked : min(time, average track duration)
%     6. Average speed in pixels/second : mean(speed(t))
%     7. Average angular speed in radians/second : mean abs(dtheta/dt (t))
%     8. Average length of object in pixels : mean (spineLength(t))
%     9. Average relative length of object, where relative length = current length / mean length :  8/mean(spineLength)
%     10. Average width of object in pixels : area / spinelength(t)
%     11. Average relative width of object : 10 / mean( area / spinelength(t))
%     12. Average aspect ratio of object (width / length) : 1/scovRatio
%     13. Average relative aspect ratio (current W/L) / (mean W/L) : 12 / mean(12)
%     14. Average “end wiggle”: angle in radians between the last 20% of the body and the rest of the
%     body (using whichever end shows a greater angle). : mean(abs(sspineTheta(t)))
%     15. Average number of pixels filled as part of object detection : mean(area(t))
% %}
% 
% data = zeros(15, length(frameNum));
% data(1,:) = round(frameNum);
% data(2,:) = et;
% data(3,:) = round(numTracked);
% data(4,:) = round(data(3,:));
% data(5,:) = min(et, mean(trackDuration));
% data(6,:) = spvstime;
% data(7,:) = dthetavstime;
% data(8,:) = slvstime;
% data(9,:) = slvstime / mean(sl(isfinite(sl)));
% data(10,:) = widthvstime;
% data(11,:) = widthvstime/mean(width(isfinite(width)));
% data(12,:) = arvstime;
% data(13,:) = arvstime/mean(1./scr(isfinite(scr)));
% data(14,:) = stvstime;
% data(15,:) = areavstime;
% if (giveUpdates)
%     disp([num2str(toc(ts1)), ' s: writing summary file']);
% end
% 
% fid = fopen([prefix '.summary'], 'wt');
% for j = frameNum
%     fprintf(fid, '%8g ', data(:,j));
%     if (any(startFrame == j) || (endFrame == j))
%         fprintf(fid, ' %%%% ');
%         if (any(startFrame == j))
%             fprintf (fid, '0 %d ', trackNum(startFrame == j));
%         end
%         if (any(endFrame == j))
%             fprintf (fid, '%d 0 ', trackNum(endFrame == j));
%         end
%     end
%     if(makeMultipleBlobsFile && any(endFrame == j))
%         fprintf (fid, ' %%%%%% ');
%         fprintf (fid, '%d %d.%d ', trackNum(endFrame == j), fileNum(endFrame == j), locInFile(endFrame == j));
%     end
%     fprintf(fid, '\n');
% end
% fclose(fid);
% if (giveUpdates)
%     disp([num2str(toc(ts1)), ' s: finished!']);
% end