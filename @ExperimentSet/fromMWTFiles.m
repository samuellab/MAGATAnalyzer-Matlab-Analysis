function eset = fromMWTFiles(basedir, camcalinfo, varargin)
% loads a set of experiments from MWT files
% eset = fromMWTFiles(basedir, camcalinfo, varargin);% 
%
% output: an Experiment Set
% inputs: basedir - a directory containing subdirectories representing MWT
%         experiments OR
%         a directory containing zip files representing MWT experiments
%camcalinfo, CameraCalibration struct or [] to ignore
%'parallel', [false]/true - whether to load in parallel using matlab
%       parallel computing library
%'interpTime', (time in seconds) - assign this value to derivation rules
% interpolation time, instead of autocomputing


existsAndDefault('camcalinfo', []);
parallel = false;
interpTime = [];
varargin = assignApplicable(varargin);
if (isempty(interpTime))
    warning('Interpolation time will be automatically determined. This may cause errors if multiple esets are later combined.');
end

eset = ExperimentSet();
if (~exist (basedir, 'file')) %checks for directories and files
    eset = repmat(eset, 0);
    disp ([basedir ' does not match any files or directories']);
    return;
end

[~,~,ext] = fileparts(basedir);
if (strcmpi (ext, '.zip') || (~isempty(dir(fullfile(basedir, '*.summary')))))
    disp ('called on mwt zip file / directory');
    fn = {basedir};
else    
    d = dir(basedir);
    d1 = d([d.isdir]);
    if (~isempty(d1))
        tf = cellfun(@(s) s(1) ~= '.', {d1.name}) & cellfun(@(s) ~isempty(dir(fullfile(basedir, s, '*.summary'))), {d1.name});
        d1 = d1(tf);
    end
    d2 = d(~[d.isdir]);
    if (~isempty(d2))
        
        tf = false(size(d2));
        for j = 1:length(d2)
            [~,~,ext] = fileparts(d2(j).name);
            tf(j) = strcmpi(ext, 'zip');
        end
        d2 = d2(tf);
    end
    d = [d1 d2];
    for j = 1:length(d)
        fn{j} = [basedir d(j).name];  %#ok<*AGROW>
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
        expt(n) = Experiment.fromMWTFile(fn{n}, camcalinfo);
        disp(['finished loading file #' num2str(n) ' : t = ' num2str(toc(ts1))]);
    end
    if (closepool)
        disp('closing parallel processes'); ts2 = tic;
        matlabpool close;
        toc(ts2);
    end
    disp(['total time to load ' num2str(length(fn)) ' files = ' num2str(toc(ts1))]);
else
    ts1 = tic;
    for n=1:length(fn)
        disp(['Loading file #' num2str(n) ' : ' fn{n}]);
        expt(n) = Experiment.fromMWTFile(fn{n}, camcalinfo);
        disp(['finished loading file #' num2str(n) ' : t = ' num2str(toc(ts1))]);
    end
end
eset.expt = expt;
try
    if (isempty(interpTime))
        et = [expt.elapsedTime];
        interpTime = percentile(diff(et), 0.1);
    end
    for j = 1:length(expt)
        expt(j).dr.interpTime = interpTime;
    end
    eset.expt.updateDerivationRules;
    %{
    mdr = num2str(interpTime);
    eset.evaluateTrackExpression(['track.dr.interpTime = ' mdr ';']);
    
    %}
    eset.executeTrackFunction('recalculateDerivedQuantities');
    eset.executeExperimentFunction('assignGlobalQuantities');
catch me
    disp(me.getReport());
end