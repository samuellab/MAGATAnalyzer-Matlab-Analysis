function eset = fromMatFiles(fstub, fileinds, segment)
% reloads experiment set from mat files
% function eset = fromMatFiles(fstub, fileinds)
%
% loads every experiment in eset from a separate mat flile
% variable name is experiment_j in file fstub_experiment_j.mat
%
% outputs: an experiment set
% inputs:
% FSTUB: everything preceding _experiment_j.mat in the file name
%      : or the path to a directory containting files
% FILEINDS: which (j) experiments to load
%   pass fileinds empty or 'all' to load all with correct form
% SEGMENT: whether to segment tracks after loading (segmentation is not
%   saved currently) - default false
% example:
% fstub = 'D:\Marc Processed\maggots\ethyl acetate 4 pct 20 2000\odor4pct'
% eset = ExperimentSet.fromMatFiles(fstub);
existsAndDefault('segment', false);

if (~existsAndDefault ('fstub', []))
    [filename,pathname] = uigetfile('*experiment*.mat', 'select mat file(s) to load','MultiSelect', 'on');
    if (~iscell(filename))
        filename = {filename};
    end
    
    %{
    fileinds = zeros(size(filename));
    temp = regexp(filename{1}, '_experiment_(\d+).mat', 'split');
    fstub = fullfile(pathname, temp{1});
    for j = 1:length(filename)
        temp = regexp(filename{j}, '_experiment_(\d+).mat', 'tokens');
        fileinds(j) = str2double(temp{1});
    end
    %}
    for j = 1:length(filename)
        filenames{j} = fullfile(pathname, filename{j});
    end
    fileinds = [];
  %  fileinds = 1:length(filenames);
else if iscell(fstub)
        filenames = fstub;
        fileinds = [];
else if isdir(fstub)
    d = dir(fullfile(fstub, '*experiment*.mat'));
    if (isempty(d))
        d = dir(fullfile(fstub, '*.mat'));
        if (isempty(d))
            disp (['no files with appropriate name structure found in ' fstub]);
            eset = [];
            return;
        end
    end
    for j = 1:length(d)
        filenames{j} = fullfile(fstub, d(j).name);
    end
    fileinds = [];
    %fileinds = 1:length(filenames);
else if (~exist('fileinds', 'var') || isempty(fileinds) || strcmpi(fileinds, 'all'))
    d = dir([fstub '_experiment_*.mat']);
    fileinds = zeros(size(d));
    for j = 1:length(d)
        ind = strfind(d(j).name, '_');
        if (~isempty(ind))
            ind = ind+1;
            ind2 = strfind(d(j).name, '.mat');
            if (~isempty(ind2))
                inds = ind(end):(ind2(end)-1);
                %d(j).name(inds)
                fileinds(j) = str2double(d(j).name(inds));
            end
        end
    end
    fileinds = fileinds(fileinds ~= 0);
end;end;end;end;

if (~exist('filenames', 'var'))
    for j = 1:length(fileinds)
        filenames{j} = [fstub '_experiment_' num2str(fileinds(j)) '.mat'];
    end
end


ts1 = tic;

eset = ExperimentSet();
foundany = false;
for k = 1:length(filenames)  
    try
        result = load (filenames{k});
        fldnm = fieldnames(result);
        nfound = 0;
        for j = 1:length(fldnm)
            if (isa(result.(fldnm{j}), 'Experiment'));  
                if (~foundany)
                    eset.expt = result.(fldnm{j});
                    foundany = true;
                else
                    eset.expt = [eset.expt result.(fldnm{j})];
                end
                nfound = nfound + 1;
            end
        end
        toc(ts1);
        if (nfound ~= 1)
            disp (['warning: loaded ' num2str(nfound) ' experiments from ' filenames{k}]);
        end
    catch me
        disp (['problem loading experiment from matfile: ' filenames{k}]);
        disp (me.getReport());
    end
end

        

if (segment)
    disp('segmenting');
    try
        eset.executeTrackFunction('segmentTrack');
    catch me
        disp(me.getReport());
        disp ('segmentation failed');
    end
end
disp (['finished - ' num2str(toc(ts1))]);
% 
% 
% if (matlabpool('size') ~= 0)
%     %fileinds
%     parfor k = 1:length(fileinds) 
%         disp(['start ' num2str(k) ' - ' num2str(toc(ts1))])
%         j = fileinds(k);
%         result = load ([fstub '_experiment_' num2str(j) '.mat'], ['experiment_' num2str(j)]);
%         fldnm = fieldnames(result);
%         expt(k) = result.(fldnm{1});
%         disp(['end ' num2str(k) ' - ' num2str(toc(ts1))]);
%     end
% 
%     eset = ExperimentSet();
%     eset.expt = expt;
%     disp(['loaded all - ' num2str(toc(ts1))]);
% else
%     for k = 1:length(fileinds) 
%         j = fileinds(k);
%         load ([fstub '_experiment_' num2str(j) '.mat'], ['experiment_' num2str(j)]);
%         toc(ts1);
%     end
% 
%     eset = ExperimentSet();
%     %we have to initialize eset.expt to be an array of experiment or we
%     %have problems later
%     j = fileinds(1);
%     eval(['eset.expt = experiment_' num2str(j) ';']);
% 
%     for k = 1:length(fileinds) 
%         j = fileinds(k);
%         eval(['eset.expt(' num2str(k) ') = experiment_' num2str(j) ';']);
% 
%     end
% end
