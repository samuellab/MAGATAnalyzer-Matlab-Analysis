function str2 = downsizeStruct (str)
%function str2 = downsizeStruct (str)
%
%converts doubles to singles to save memory space


fn = fieldnames (str);
for k = 1:length(str)
    for j = 1:length(fn)
        if (isa(str(k).(fn{j}), 'float'))
            str2(k).(fn{j}) = single(str(k).(fn{j}));
        else
            if isa(str(k).(fn{j}), 'struct')
                str2(k).(fn{j}) = downsizeStruct(str(k).(fn{j}));
            else
                str2(k).(fn{j}) = str(k).(fn{j});
            end
        end
    end
end