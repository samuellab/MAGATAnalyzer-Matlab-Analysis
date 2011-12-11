function data = toHashMap(obj, varargin)
%convert a generic object or structure into a java hashmap
%function data = toHashMap(obj, varargin)
%
%converts an object or structure into a java hashmap, where
%each field name is a key, and each field value is a value
%if the field value is itself an object, then that object is converted to
%a hashmap, and that hashmap is stored as the value
%
%to prevent infinite recursion, if the a field is a handle object, it is
%only converted to a hashmap if handleRecurseLevel > 0 (default 0)
%in this case, handleRecurseLevel is decremented before being passed to 
%toHashMap(subfield)
%
%if obj is an array of objects, a java array of hashmaps is returned
handleRecurseLevel = 0;
varargin = assignApplicable(varargin);
data = java.util.HashMap();
if (length(obj) > 1)
    arr = javaArray('java.util.HashMap', length(obj));
    for j = 1:length(obj)
        arr(j) = toHashMap(obj(j), 'handleRecurseLevel', handleRecurseLevel, varargin{:});
    end
    data = arr;
    return;
end

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
        data.put(f{j}, toHashMap(obj.(f{j}), 'handleRecurseLevel', handleRecurseLevel, varargin{:}));
    else
        data.put(f{j}, matlabArrToJavaArr(obj.(f{j})));
    end
    
end

function javatype = getJavaType (matlabarr)
matlabtypes = {'logical', 'single', 'double', 'char', 'int8', 'uint8', 'int16', 'uint16', 'int32', 'uint32', 'int64', 'uint64'};
javatypes = {'Boolean', 'Float', 'Double', 'String', 'Byte', 'Byte', 'Short', 'Short', 'Integer', 'Integer', 'Long', 'Long'};
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
%            warning ('could not convert cell to array');
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
matlabarr = matlabarr(:);
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