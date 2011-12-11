function eset = fromFiles (varargin)
% loads a set of experiments from .bin files
% function eset = ExperimentSet.fromFiles (varargin)
% 
% Load a bunch of tracks.bin into the experiment. modified from ARC code
%
% output: an Experiment Set
% inputs: variable
% 'minpts', XX sets the minimum number of points in the tracks to load
% 'loadcontour', true/false (default true); whether to load Maggot ctrs
% with no (additional) arguments, user selects bin files to load
% if 1st argument is a directory, loads all .bin files from that
% directory; ignores all other args
% with >1 argument, each argument is the file name to load
%'camcalinfo', CameraCalibration struct or [] to ignore
%'parallel', [false]/true - whether to load in parallel using matlab
%       parallel computing library
%'sortbydate', [false]/true - if expts should be sorted by file date;  if
%false they are the default sort order, which I think is alphabetical

minpts = 100;
loadcontour = true;
camcalinfo = [];
parallel = false;
sortbydate = false;
varargin = assignApplicable(varargin);
eset = ExperimentSet();
if (isempty(varargin))
    [fn, basedir] = uigetfile('*.bin','Select the .bin files you would like to use', 'MultiSelect', 'on');
    if (~iscellstr(fn))
        fn = {fn};
    end
    for j = 1:length(fn)
        fn{j} = [basedir fn{j}];
    end
else
    if (isdir(varargin{1}))
        basedir = varargin{1};
        if (basedir(end) ~= '\')
            basedir = [basedir '\'];
        end
        d = dir([basedir '*.bin']);
        if (sortbydate)
            [~,I] = sort([d.datenum]);
            d = d(I);
        end
        for j = 1:length(d)
            fn{j} = [basedir d(j).name];  %#ok<*AGROW>
        end
    else
        fn = varargin;
    end
end


if (parallel && length(fn) > 1)
    ts1 = tic;
    if (matlabpool('size') == 0)
        disp ('opening parallel processes'); 
        matlabpool;
        toc(ts1);
        closepool = true;
    else
        closepool = false;
    end
    
    
    parfor n=1:length(fn)
        disp(['Loading file #' num2str(n) ' : ' fn{n}]);
        ind = strfind(fn{n}, '.');
        if (isempty(ind))
            disp(['problem with ' fn{n}]);
        end
        %timfn=[fn{n}(1:ind(end)) 'tim'];
        if (length(camcalinfo) >= n)
            expt(n) = Experiment.fromFile (fn{n}, [], loadcontour, camcalinfo(n), minpts);
        else
            expt(n) = Experiment.fromFile (fn{n}, [], loadcontour, camcalinfo, minpts);
        end
        disp(['finished loading file #' num2str(n) ' : t = ' num2str(toc(ts1))]);
    end
    if (closepool)
        disp('closing parallel processes'); ts2 = tic;
        matlabpool close;
        toc(ts2);
    end
    disp(['total time to load ' num2str(length(fn)) ' files = ' num2str(toc(ts1))]);
else
    for n=1:length(fn)
        disp(['Loading file #' num2str(n) ' : ' fn{n}]);
        ind = strfind(fn{n}, '.');
        if (isempty(ind))
            disp(['problem with ' fn{n}]);
        end
        %timfn=[fn{n}(1:ind(end)) 'tim'];
         if (length(camcalinfo) >= n)
            expt(n) = Experiment.fromFile (fn{n}, [], loadcontour, camcalinfo(n), minpts);
        else
            expt(n) = Experiment.fromFile (fn{n}, [], loadcontour, camcalinfo, minpts);
        end
    end
end
eset.expt = expt;
