function n = objNargout (obj, func)
%function n = objNargout (obj, func)

%get metaclass information
if (ischar(obj))
    mc = eval(['? ' obj]);
else
    mc = metaclass(obj);
end

%get all methods
meth = [mc.Methods{:}];

ind = find(strcmp(func, {meth.Name}), 1);

if (isempty(ind))
    n = [];
    return;
end

if (isempty(meth(ind).OutputNames))
    n = 0;
    return;
end

if any(strcmp('varargout', meth(ind).OutputNames))
    n = -1;
    return;
end

n = length(meth(ind).OutputNames);