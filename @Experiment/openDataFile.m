function openDataFile(expt)
%opens the data file associated with expt for reading
%function openDataFile(expt)
%
%opens the data file associated with expt for reading (e.g. for reloading
%images)
%inputs: EXPT, a member of the Experiment class
if (expt.fid ~= 0)
    try
        fclose(expt.fid);
    catch me
        % do nothing, it was probably already closed
    end
end

expt.fid = fopen(expt.fname, 'r');