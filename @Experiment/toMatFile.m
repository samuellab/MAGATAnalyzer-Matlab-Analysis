function toMatFile(experiment, varargin)
%function toMatFile(experiment, varargin)
%
%saves an experiment to a .mat file with the same file name as the .bin
%file, but with a different extension
%extradir, 'extra\directories' -- places matfile in a subdirectory 
extradir = '';
tracksInSubDir = false;
fname = '';
varargin = assignApplicable(varargin);
if (isempty(fname))
    [p, f] = fileparts(experiment.fname);
    if (~isempty(extradir))
        fname = fullfile(p, extradir, [f '.mat']);
    %    save(fullfile(p, extradir, [f '.mat']), 'experiment')
    else
        fname = fullfile(p, [f '.mat']);
    %    save(fullfile(p, [f '.mat']), 'experiment');
    end
end
if (tracksInSubDir)
    experiment = Experiment(experiment); %clone to avoid sending changes backwards
    [p,f] = fileparts(fname);
    experiment.savedTrackRelDir = [shortenFileStub(f) ' - tracks'];
    experiment.savedTrackDir = fullfile(p, experiment.savedTrackRelDir);
    if (~exist(experiment.savedTrackDir, 'dir'))
        mkdir(experiment.savedTrackDir);
    end
    for j = 1:length(experiment.track)
        experiment.track(j).toMatFile(fullfile(experiment.savedTrackDir, ['track' num2str(j) '.mat']));
    end
    experiment.track = [];
end

save(fixFileNameWin(fname), 'experiment');