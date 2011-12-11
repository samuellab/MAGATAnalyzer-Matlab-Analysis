function str = objToYaml(key, obj, handleRecurseLevel)
%parses an object to yaml using snakeyaml
%function str = objToYaml(obj)
%
existsAndDefault('handleRecurseLevel', 0);

import('org.yaml.snakeyaml.Yaml');
data = java.util.HashMap();
arr = javaArray('java.util.HashMap', length(obj));
for j = 1:length(obj)
    arr(j) = objToHashmap(obj(j), handleRecurseLevel);
end
data.put(key, arr);
yaml = Yaml();
str = yaml.dump(data);

%{

function data = objToHashmap(obj, handleRecurseLevel)
data = java.util.HashMap();
if (length(obj) > 1)
    arr = javaArray('java.util.HashMap', length(obj));
    for j = 1:length(obj)
        arr(j) = objToHashmap(obj(j), handleRecurseLevel);
    end
    data = arr;
    return;
end

existsAndDefault('handleRecurseLevel', 0);
if (isa(obj, 'handle'))
    if (handleRecurseLevel <= 0)
        warning ('oth:mrl', 'maximum handle recursion reached');
        data.put('handle type', class(obj));
        return;
    else
        handleRecurseLevel = handleRecurseLevel - 1;
    end
end
f = fieldnames(obj);
for j = 1:length(f)
    if (isobject(obj.(f{j})) || isstruct(obj.(f{j})))
        data.put(f{j}, objToHashmap(obj.(f{j}), handleRecurseLevel));
    else
        data.put(f{j}, matlabArrToJavaArr(obj.(f{j})));
    end
    
end

function javatype = getJavaType (matlabarr)
matlabtypes = {'logical', 'single', 'double', 'char', 'int8', 'uint8', 'int16', 'uint16', 'int32', 'uint32', 'int64', 'uint64'};
javatypes = {'Boolean', 'Float', 'Double', 'String', 'Byte', 'Byte', 'Short', 'Short', 'Int', 'Int', 'Long', 'Long'};
for j = 1:length(matlabtypes)
    if (isa (matlabarr, matlabtypes{j}))
        javatype = ['java.lang.' javatypes{j}];
        return;
    end
end
if (isnumeric(matlabarr))
    javatype = 'java.lang.Double';
    return;
end
javatype = 'java.lang.Object';

function arr = matlabArrToJavaArr (matlabarr)
if ischar(matlabarr)
    matlabarr = {matlabarr};
end
if (iscell(matlabarr))
    if (iscellstr(matlabarr))
        javatype = 'java.lang.String';
        arr = javaArray('java.lang.String', max(1,length(matlabarr)));
    else
        try 
            matlabarr = cell2mat(matlabarr);
            javatype = getJavaType(matlabarr);
        catch
            warning ('could not convert cell to array');
            javatype = 'java.lang.Object';
        end
    end
else
    javatype = getJavaType(matlabarr);
end
arr = javaArray(javatype,  max(1,length(matlabarr)));

if (iscellstr(matlabarr))
    for j = 1:length(matlabarr)
        arr(j) = javaObject(javatype, matlabarr{j});
    end
    return;
end
if (iscell(matlabarr))
    for j = 1:length(matlabarr)
        if isa(matlabarr{j}, 'function_handle')
            msg = func2str(matlabarr{j});
            arr(j) = java.lang.String(msg);
        else
            arr(j) = matlabArrToJavaArr(matlabarr{j});
        end
    end
    return;
end

for j = 1:length(matlabarr)
    if (strcmpi(javatype, 'java.lang.Object'))
        if isa(matlabarr(j), 'function_handle')
            msg = func2str(matlabarr(j));
        else
            msg = 'unknown type';
        end
        arr(j) = java.lang.String(msg);
    else
        arr(j) = javaObject(javatype,matlabarr(j));
    end
end
%}