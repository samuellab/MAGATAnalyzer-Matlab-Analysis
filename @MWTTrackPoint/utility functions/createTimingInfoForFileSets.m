function createTimingInfoForFileSets (batfilename, varargin)
%function createTimingInfoForFileSets (batfilename, varargin)
%
%looks through a bat file, detects input & output files, then adds timing
%info files (.tim files) if these do not already exist
%
%optional param, value pairs
%'verbose', true/false -- whether or not to display additional info as
%processing proceeds
%'extension', '.ext' -- extension if not jpg, see createTimFile

verbose = false;
varargin = assignApplicable(varargin);

if (verbose)
    ts = tic;
    disp (['reading in ' batfilename]);
end

str = fileread(batfilename);
[directories, fstubs, outnames] = getFileInfoFromBatFile(str);
if (length(directories) ~= length(fstubs) || length(fstubs) ~= length(outnames))
    disp ('Error! Did not parse correctly')
    directories
    fstubs
    outnames
end

if (verbose)
    toc (ts);
    disp (['detected ' num2str(length(directories)) ' file sets']);
end

for j = 1:length(outnames)
    ind = find(outnames{j} == '.',1,'last');
    outnames{j} = [outnames{j}(1:ind) 'tim'];
end

for j = 1:length(outnames)
    if (exist(outnames{j}, 'file'))
        if (verbose)
            disp ([outnames{j} ' already exists; taking no action']);
        end
        continue;
    else
        createTimFile(directories{j}, fstubs{j}, outnames{j}, 'verbose', verbose, varargin{:});
    end
end
