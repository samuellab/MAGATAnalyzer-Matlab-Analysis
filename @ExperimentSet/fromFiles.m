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

%Set params
minpts = 100;
loadcontour = true;
camcalinfo = [];
parallel = false;
sortbydate = false;
fixedInterpTime = [];
varargin = assignApplicable(varargin);

%Build the eset
eset = ExperimentSet();

%Get the file names
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

%Ignore any files that have been marked as bad
valid = true(size(fn));
for j = 1:length(fn)
    [bd,f]  = fileparts(fn{j});
    if exist (fullfile(bd, [f '.bad']), 'file')
        warning ('EFF:BAD', ['file: ' f ' is marked as bad and will not be loaded. To load file, remove ' f '.bad from directory and rerun command']);
        valid(j) = false;
    end
end
fn = fn(valid);
valid = true(size(fn));


expt = repmat(Experiment(), size(valid));

if (parallel && length(fn) > 1)
%Load the experiments in parallel
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
        try 
            disp(['Loading file #' num2str(n) ' : ' fn{n}]);
            ind = strfind(fn{n}, '.');
            if (isempty(ind))
                disp(['problem with ' fn{n}]);
            end
            %timfn=[fn{n}(1:ind(end)) 'tim'];
            if (length(camcalinfo) >= n)
                %Actually load the experiment!!
                %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                [temp,valid(n)] = Experiment.fromFile (fn{n}, [], loadcontour, camcalinfo(n), minpts);
                %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                if (valid(n) && ~isempty(temp))
                    expt(n) = temp;
                end
            else
                %Actually load the experiment!!
                %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                [temp,valid(n)] =  Experiment.fromFile (fn{n}, [], loadcontour, camcalinfo, minpts);
                %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                if (valid(n)&& ~isempty(temp))
                    expt(n) = temp;
                end
            end
            disp(['finished loading file #' num2str(n) ' : t = ' num2str(toc(ts1))]);
        catch me
            %valid(n)=0;
            disp('**********');
            disp(me.getReport());
            disp([fn{n} 'failed to process!']);
            disp('**********');
        end
    end
    if (closepool)
        disp('closing parallel processes'); ts2 = tic;
        matlabpool close;
        toc(ts2);
    end
    disp(['total time to load ' num2str(length(fn)) ' files = ' num2str(toc(ts1))]);
    
else
%Load the experiments in series:
    for n=1:length(fn)
        disp(['Loading file #' num2str(n) ' : ' fn{n}]);
        ind = strfind(fn{n}, '.');
        if (isempty(ind))
            disp(['problem with ' fn{n}]);
        end
        %timfn=[fn{n}(1:ind(end)) 'tim'];
         if (length(camcalinfo) >= n)
             %Use separate camcalinfo's for each expt
             try
                [ppp, fff] = fileparts(fn{n});
                if (exist(fullfile(ppp, [fff '.bad']), 'file'))
                    disp ([fn{n} ' is marked as a bad file -- ignoring.']);
                    disp (['fix problem and delete file: ' fff '.bad']);
                    temp = [];
                    valid(n) = false;
                else
                    %Actually load the experiment!!
                    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                    temp = Experiment.fromFile (fn{n}, [], loadcontour, camcalinfo(n), minpts);
                    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                end
                if (~isempty(temp))
                    expt(n) = temp;
                else
                    valid(n) = false;
                end
             catch me
                 disp('**********');
                 disp([fn{n} 'failed to process!']);
                 disp('**********');
                 disp (me.getReport());
                 valid(n) = false;
                 [ppp, fff] = fileparts(fn{n});
                 fidd = fopen(fullfile(ppp, [fff '.bad']),'wt');
                 fprintf(fidd, me.getReport());
                 fclose(fidd);
             end
         else
            %Use the same camcalinfo for all expt's
            try
                %Actually load the experiment!!
                %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                temp = Experiment.fromFile (fn{n}, [], loadcontour, camcalinfo, minpts);
                %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                if (~isempty(temp))
                    expt(n) = temp;
                else
                    valid(n) = false;
                end
             catch me
                 disp('**********');
                 disp([fn{n} 'failed to process!']);
                 disp('**********');
                 disp (me.getReport());
                 valid(n) = false;
             end
        end
    end
end

if (~exist ('expt', 'var') || isempty(expt))
    disp ('no valid experiments loaded');
    eset = repmat(eset, 0);
    return;
end

%Only keep the valid experiments
eset.expt = expt(valid);
try %restore try deleted by nbernat 05-20 -- don't understand how this worked at all afterwards
    %Calculate derived quantities
    dr = [expt.dr];
    if ((~isempty(fixedInterpTime) && fixedInterpTime > 0) || min([dr.interpTime]) < 0.98 * max([dr.interpTime])) %#ok<BDSCI>
        if ((~isempty(fixedInterpTime) && fixedInterpTime > 0))  %#ok<BDSCI>
            mdr = num2str(fixedInterpTime);
        else
            mdr = num2str(min([dr.interpTime]));
        end
        eset.evaluateTrackExpression(['track.dr.interpTime = ' mdr ';']);
        mdr = min([dr.interpTime]);
        for j = 1:length(expt)
            expt(j).dr.interpTime = mdr;
        end
        eset.executeTrackFunction('recalculateDerivedQuantities');
        eset.executeExperimentFunction('assignGlobalQuantities');
    end
catch me
    disp(me.getReport());
end