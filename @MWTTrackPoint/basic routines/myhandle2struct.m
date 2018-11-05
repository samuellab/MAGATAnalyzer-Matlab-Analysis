function s = handle2struct(obj, parents)
%function s = handle2struct(obj, parents)
%
%converts all handles to a structure in order to get at true memory size
%to avoid infinite recursion, we won't call handle to struct on any handle 
%that has already been called
existsAndDefault('parents', repmat(handle, [0 0]));
try
    warning('off', 'all');
    if (isa(obj, 'handle'))
        if (any(parents == obj))
            s = obj;
            return;
        end
        parents = [parents, obj];
        s = struct(obj);
    else
        s = obj;
    end
catch me
    me.getReport();
    obj
    parents
end
if (isa (s, 'struct'))
    f = fieldnames(s);
    for j = 1:length(f)
        for k = 1:length(s.(f{j}))
            s.(f{j})(k) = handle2struct(s.(f{j})(k),parents);
        end
    end
end
