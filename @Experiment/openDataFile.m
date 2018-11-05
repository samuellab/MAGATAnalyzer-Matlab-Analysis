function openDataFile(expt)
%opens the data file associated with expt for reading
%function openDataFile(expt)
%
%opens the data file associated with expt for reading (e.g. for reloading
%images)
%inputs: EXPT, a member of the Experiment class
if (isa(expt.track(1), 'MWTTrack'))
    %no data file to read in same way for multi worm tracker
    return;
end

if (expt.fid ~= 0)
    try
        fclose(expt.fid);
    catch me
        % do nothing, it was probably already closed
    end
end

expt.fid = fopen(expt.fname, 'r');