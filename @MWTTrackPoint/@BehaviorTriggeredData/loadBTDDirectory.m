function btdstruct = loadBTDDirectory (basedir, varargin)
%function btdstruct = loadBTDDirectory (basedir, varargin)


if (iscell(basedir))
    if (length(basedir) == 1)
        btdstruct = BehaviorTriggeredData.loadBTDDirectory (basedir{1}, varargin);
        return;
    end
    for j = 1:length(basedir)
        bs(j) = BehaviorTriggeredData.loadBTDDirectory (basedir{j}, varargin{:});
    end
    btdstruct.srcdir = {bs.srcdir};
    btdstruct.btd = [bs.btd];
    return;
end

fileList = {};
varargin = assignApplicable(varargin);

if (isempty(fileList))
    d = dir(fullfile(basedir, 'btd_*.mat'));
    f = {d(~[d.isdir]).name};
else
    f = fileList;
    if ~iscell(f)
        f = {f};
    end
end
for j = 1:length(f)
   s = load(fullfile(basedir, f{j}), 'btd');
    btd(j) = s.btd; %#ok<AGROW>
    valid(j) = ~isempty(btd(j).glt); %#ok<AGROW>
end
if (isempty(f))
    error (['bad directory ' basedir]);
end
if (~all(valid))
    inds = find(~valid);
    for j = 1:length(inds)
        disp ([f{inds(j)} ' had problems loading']);
    end
end


btd = btd(valid);
btdstruct.srcdir = basedir;
btdstruct.btd = btd;




