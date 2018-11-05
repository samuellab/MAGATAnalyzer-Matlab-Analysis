function [ptType, valid] = getPtType (fname)
%loads an experiment from a bin file
%[expt, valid] = fromFile (fname, timfname, loadContour, camcalinfo, minTrackLength)
%this is a static method of the Experiment class (Experiment.fromFile)
%
%outputs: 
%EXPT, a member of the experiment class
%VALID, whether the experiment loaded without errors
%inputs:
%FNAME: name of .bin file to load
%TIMFNAME: timining information file (.tim);
%   default: change extension of fname to .tim
%LOADCONTOUR: whether to load the contour if this is a maggot track:
%   default TRUE
%CAMCALINFO: camera calibration struct (ask Marc); pass empty ([]) to ignore
%   default: []
%MINTRACKLENGTH: minimum length of a track (in points) to load from disk
%   default: 1
tic
try 
        
    valid = true;
    
    
    
    fid = fopen(fname, 'r');
    code = fread(fid, 1, 'int32');
    fclose(fid);
    
    switch (bitshift(code, -16))
        case 1
            ptType = TrackPoint();
        case 2
            ptType = ImTrackPoint();
        case 3
            ptType = OldMaggotTrackPoint();
            
        case 4
            ptType = MaggotTrackPoint();
            
        otherwise
            disp('invalid code: I don''t know what kind of point I''m loading');
            ptType = [];
            valid = false;
            
    end
catch me
    
    disp(me.getReport());
    valid = false;
    ptType = [];
    if (exist('fid', 'var'))
        fclose(fid);
    end
end
    