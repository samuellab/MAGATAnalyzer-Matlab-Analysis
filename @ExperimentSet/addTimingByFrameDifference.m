function addTimingByFrameDifference(eset, deltaT)
%function addTimingByFrameDifference(eset, deltaT)
%
%iterates through experiments and assigns timing to any that do not already
%have elapsedTime defined

for j = 1:length(eset.expt)
    if (isempty(eset.expt(j).elapsedTime)) 
        try
            [~,nm] = fileparts(eset.expt(j).fname);
            warning('ESET:NOTIME', [nm ': no timing information -- adding guessed info']);
        catch
        end
        eset.expt(j).setTimingByFrameDifference(deltaT, false);
    end
end