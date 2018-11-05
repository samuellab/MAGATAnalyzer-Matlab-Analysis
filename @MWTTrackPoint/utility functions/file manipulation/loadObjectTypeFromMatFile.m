function objs = loadObjectTypeFromMatFile (fname, objType)
%function loadObjectTypeFromMatFile (fname, objType)
%
%only loads objects of a specified type from file
%if multiple objects are present in the file, attempts to concatenate them
%as an array
%if that fails, returns a cell array
if (nargin < 2)
    error ('need object type');
end
try
    s = load(fixFileNameWin(fname));
catch me
    disp ('failed to read from file');
    disp (me.getReport());
    objs = [];
    return;
end

fn = fieldnames(s);
for j = 1:length(fn)
    if ~isa(s.(fn{j}), objType)
        s = rmField(s, fn{j});
    end
end

fn = fieldnames(s);
for j = 1:length(fn)
    objs{j} = s.(fn{j}); %#ok<AGROW>
end

try
    objs = [objs{:}];
catch
    %intentionally blank
end