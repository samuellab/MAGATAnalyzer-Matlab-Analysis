function tr = detectTimingGaps (basedir, varargin)
% looks through all .tim files and detects large gaps in acquistition time
% function timingreport = detectTimingGaps (basedir, varargin)
% optional arguments:
%       mingap = 1000; minimum gap time in milliseconds

mingap = 1000;
varargin = assignApplicable(varargin);

timfiles = recursiveDirectorySearch(basedir, '*.tim');
if (isempty(timfiles))
    tr = [];
    return;
end

for j = 1:length(timfiles)
    data = load(timfiles{j});
    t = data(:,1);
    ind = find(diff(t) < 0, 1, 'first');
    while (~isempty(ind))
        t((ind+1):end) = t((ind+1):end) + 65536;
        ind = find(diff(t) < 0, 1, 'first');
    end
    tr(j).fname = timfiles{j};
    tr(j).maxdt = max(diff(t));
    tr(j).nbad = sum(diff(t) > mingap);
    tr(j).haserrors = tr(j).nbad > 0;
    tr(j).et = t;
end

disp ([num2str(length(timfiles)) ' files examined. ' num2str(sum([tr.haserrors])) ' have large timing gaps.']);
if (any([tr.haserrors]))
    disp ('the following files have errors: ');
    for j = 1:length(tr)
        if (tr(j).haserrors)
            disp([tr(j).fname ' - max gap ' num2str(tr(j).maxdt)]);
        end
    end
end