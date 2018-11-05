function mystruct = stripFunctionHandles( mystruct )
% function mystruct = stripFunctionHandles( mystruct )
%   converts any function handles in struct or nested structures
%   to strings
%   the purpose is to remove stack data saved with the function that leads
%   to large file sizes

if (isstruct(mystruct))
    fn = fieldnames(mystruct);
    for k = 1:length(mystruct)
        for j = 1:length(fn)
            mystruct(k).(fn{j}) = stripFunctionHandles(mystruct(k).(fn{j}));
        end
    end
    return;
end
if (iscell(mystruct))
    for k = 1:length(mystruct)
        mystruct{k} = stripFunctionHandles(mystruct{k});
        
    end
    return;
end

if isa(mystruct, 'function_handle')
    mystruct = func2str(mystruct);
end

end

