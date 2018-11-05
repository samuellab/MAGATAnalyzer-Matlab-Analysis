function [expt, valid] = fromJava(jEx, fname, timfname, loadContour, camcalinfo, minTrackLength)
% Replacement for Experiment.fromFile;
%   Opens [JavaTrackExtraction.Experiment jEx]<-FILENAME
%   Converts jEx to [Experiment ex] that is consistent with the experiments
%       loaded with Ex.fromFile

% NOTE: See C:\Users\Natalie\Documents\GitHub\Matlab-Track-Analysis\@Experiment\fromFile.m
%       for processing messages

debug = true;

try 
    
    % Set flags/parameters for later processing
    if (~exist ('loadContour', 'var') || isempty ('loadContour'))
        loadContour = true;
    end
    if (~exist ('camcalinfo', 'var'))   
        camcalinfo = [];
    end
    if (~exist ('minTrackLength', 'var') || isempty(minTrackLength))
        minTrackLength = 1;
    end
    
%     disp('~~~~~~~~~ Experiment.fromJava ~~~~~~~~~');
    
    % Create the experiment
%     disp('**************************');
%     disp(['Processing jav->mat: ' fname]);
%     disp('**************************');
    expt = Experiment();
    valid = true;
    % set file name
    expt.fname = fname;
    % set camcalinfo
    expt.camcalinfo = camcalinfo;

    %...Not checking whether or not fname is a MWTTrack

    % Determine number of tracks nTracks
    %%%nTracks = jEx.getNumTracks;
    nTracks = javaMethod('getNumTracks','TrackExtractionJava.Experiment',fname);
    
    % Determine point type & allocate track array
    % ~~ 3 LTP(=BTP) ~~ 2 MTP ~~ 1 ITP ~~ 0 TP
    %%%typeCode = jEx.getTypeCode;
    typeCode = javaMethod('getPointType','TrackExtractionJava.Experiment',fname);
    switch(typeCode)
        case 3
            ptType = LarvaTrackPoint();
            tracks(nTracks) = LarvaTrack();
        case 2
            ptType = MaggotTrackPoint();
            tracks(nTracks) = MaggotTrack();
        case 1
            ptType = ImTrackPoint();
            tracks(nTracks) = Track();
        case 0
            ptType = TrackPoint;
            tracks(nTracks) = Track();
        otherwise
            disp('Invalid point type code');
            valid = false;
            [ppp, fff] = fileparts(fname);
            fidd = fopen(fullfile(ppp, [fff '.bad']),'wt');
            fprintf(fidd, 'invalid code: I don''t know what kind of point I''m loading\n');
            fclose(fidd);
            return;
    end

    
%     nPts = 0;
%     for i=0:(nTracks-1)
%         t = javaMethod('getTrack', 'TrackExtractionJava.Experiment', i,fname);
%         nPts = nPts + t.getNumPoints();
%     end
    % Loop over tracks
    disp(['Loading ' int2str(nTracks) ' tracks...']);% (' int2str(nPts) ' pts)
    ts = tic; 
    ptsLastLoaded = 0;
    lastelapsed = 0;
    reportEvery = 60;
    for i = 0:(nTracks-1)
        % ~~ Display processing messages
        elapsed = toc(ts);
        if (elapsed - lastelapsed > reportEvery)
            ptsPerMin = (60.0*ptsLastLoaded)/(elapsed - lastelapsed);
            lastelapsed = elapsed;
            disp ([num2str(elapsed) 's: ' num2str(i) '/' num2str(nTracks) ' tracks (' int2str(ptsPerMin) ' pts/min)']);
            ptsLastLoaded=0;
        end
        % ~~ Add track TrackFromJava(jTrack(i))
        disp(['track ' num2str(i)]);
        jTr = javaMethod('getTrack', 'TrackExtractionJava.Experiment', i,fname);
        tracks(i+1) = Track.fromJava(jTr, ptType, i, [], loadContour, camcalinfo, minTrackLength);
        ptsLastLoaded = ptsLastLoaded + jTr.getNumPoints();
        %%%tracks(i+1) = Track.fromJava(jEx.getTrackFromInd(int32(i)), ptType, i, [], loadContour, camcalinfo, minTrackLength);
    end
    disp(['...done loading tracks (' int2str(toc(ts)) 'sec)']);
    
    % Show a report (or don't)
    showReport = false;
    if (showReport)
        semilogy (0:10:max([tracks.npts]), hist([tracks.npts], 0:10:max([tracks.npts])), 'b-'); hold on;
        title ('distribution of track length in frames');
    end
    
    disp('Tidying & bookkeeping...');
    
    % Remove short tracks
    expt.track = tracks([tracks.npts] > minTrackLength);
    % Again, show a report (or don't)
    if (showReport)
        semilogy (0:10:max([expt.track.npts]), hist([expt.track.npts], 0:10:max([expt.track.npts])), 'r-'); hold off;
        legend ('pre trim', 'post trim');
    end
    
    % Assign ex to all track.ex's
    [expt.track.expt] = deal(expt);
    
    disp('...done');
    
    if (debug)
        try
            expt.toMatFile('extradir', 'noTime');
        catch
            disp('error saving expt (pre-timing info)')
        end
    end
        
    % Add timing info
    try
        disp ('Adding timing info...')
        expt.addtime(timfname);
        %expt.addtime(timfname);
        disp('...done adding timing info')
    catch e
        disp(e.getReport());
        disp('...done adding timing info')
    end

catch me
    valid = false;
    [ppp, fff] = fileparts(fname);
    fidd = fopen(fullfile(ppp, [fff '.bad']),'wt');
    fprintf(fidd, me.getReport());
    fclose(fidd);
    rethrow(me);
end
% disp('**************************');
% disp(['Done processing jav->mat: ' fname]);
% disp('**************************');
end