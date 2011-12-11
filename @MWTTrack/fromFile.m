function  mt = fromFile(fname, camcalinfo) 
% loads a MWT Track from file, specified by fname 
% mt = fromFile(fname)
% if fname ends with "blobs," multiple tracks are loaded
existsAndDefault('camcalinfo', []);
if (iscell(fname))
    mt = MWTTrack.fromFile(fname{1}, camcalinfo);
    for j = 2:length(fname)
        mt = [mt MWTTrack.fromFile(fname{j}, camcalinfo)];
    end
    disp ('fixing head tail orientation'); ts1 = tic;
    for j = 1:length(mt)
        try
            mt(j).fixHTOrientation;
        catch me
            j
            disp(me.getReport);
        end
            
    end
    toc(ts1);
    return;
end
[~,~,ext] = fileparts(fname);

if (strcmpi(ext, '.blob'))
    fid = fopen(fname,'rt');
    if (fid < 0 || isempty(fid))
        error (['could not open file: ' fname]);
    end
    mt = fromBlobFile(fid, camcalinfo);
    mt.fname = fname;
    fclose(fid);
    return
end
if (strcmpi(ext, '.blobs'))
    fid = fopen(fname, 'rt');
    if (fid < 0 || isempty(fid))
        error (['could not open file: ' fname]);
    end
    mt = fromBlobsFile(fid, camcalinfo);
    [mt.fname] = deal(fname);
    fclose(fid);
    return
end

warning ('Multi Worm Tracker tracks must be loaded from .blob or .blobs files');
mt = [];
return


end


function mt = fromBlobFile(fid,camcalinfo)
    mt = MWTTrack;
    mt.locInFile = fpos;
    pt = readTrackPoints(fid);
    mt.pt = pt;
    mt.startFrame = pt(1).ind;
    mt.endFrame = pt(end).ind;
    mt.npts = length(pt);    
end

function mt = fromBlobsFile(fid,camcalinfo)
    nt = 0;
    fpos = ftell(fid);
    fseek(fid, 0, 'eof');
    totalsize = ftell(fid);
    fseek(fid, fpos, 'bof');
    ts = tic;
    lastelapsed = 0;
    reportEvery = 30;
    
    while (~feof(fid) && isempty(ferror(fid)))
        fpos = ftell(fid);
        elapsed = toc(ts);
        if (elapsed - lastelapsed > reportEvery)
            lastelapsed = elapsed;
            disp ([num2str(elapsed) 's: ' num2str(ftell(fid)) '/' num2str(totalsize) ' bytes (' num2str(100*ftell(fid)/totalsize, 2) '%) loaded' ...
                num2str(elapsed*(totalsize - ftell(fid))/ftell(fid)) ' s remain']);
        end

        nt = nt+1;
        mt(nt) = MWTTrack; %#ok<*AGROW>
        mt(nt).locInFile = fpos;
        str = fgetl(fid);
        if (str(1) ~= '%' && ~feof(fid) && ~ferror(fid) && str(1) ~= -1)
            warning (['expected % to begin track designation at ' fpos ' into file']);
            disp ('discarding the following: ')
            while (str(1) ~= '%')
                disp(str);
                str = fgetl(fid);
            end
        end
        mt(nt).trackNum = sscanf(str(2:end),'%d');
        pt = readTrackPoints(fid,camcalinfo);
        mt(nt).pt = pt;
        if (isempty(pt))
            valid(nt) = false;
            continue;
        end
        valid(nt) = true;
        mt(nt).startFrame = pt(1).ind;
        mt(nt).endFrame = pt(end).ind;
        mt(nt).npts = length(pt);    
        
    end
    mt = mt(valid);
end

function mtp = readTrackPoints(fid,camcalinfo)
    fpos = ftell(fid);    
    str = fgetl(fid);
    npts = 0;
    mwtp = MWTTrackPoint;
  %  mtp = repmat(mwtp, [1 1E4]);
   % valid = false(size(mtp));
    while (~feof(fid) && isempty(ferror(fid)) && str(1) ~= -1 && str(1) ~= '%')
        try
            npts = npts+1;
            pt = mwtp.fromMWTString(str,camcalinfo);
            pt.locInFile = fpos;
            mtp(npts) = pt;
            valid(npts) = true;
        catch
            mtp(npts) = MWTTrackPoint;
            valid(npts) = false;
        end
        fpos = ftell(fid);
        str = fgetl(fid);
    end
    if (str(1) == '%')
        fseek(fid, fpos,'bof');
    end
    mtp = mtp(valid);
end

        