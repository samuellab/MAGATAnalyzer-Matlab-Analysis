function obj = objFromHashmap(hm)
%convert a java hashmap into a structure
%function obj = objFromHashmap(hm)
%
%converts a java hashmap into a class or stucture given by obj, where
%each field name is a key, and each field value is a value
%if the field value is itself an object, then that object is converted to
%a hashmap, and that hashmap is stored as the value
%
%in general, one should not expect that objFromHashmap(objToHashmap(obj))
%is the same as obj
%
%if hm is an array of hashmaps (hm is a java.util.ArrayList), obj will be
%an array of objects

if isa(hm, 'java.util.ArrayList')
    it = hm.iterator;
    for j = 1:hm.size()
        try
            no = it.next; %matlab autoconverts basic types
            if (isjava(no))                
                obj(j) = objFromHashmap(no);
            else
                obj{j} = no;
            end
        catch me
            disp('error');
            me.getReport
            break;
        end
    end
    if (iscell(obj) )
        if (length(obj) == 1)
            obj = obj{1};
        else
            if (~iscellstr(obj))
                try 
                    obj = cell2mat(obj);
                catch
                    %nothing
                end
            end
        end
    end 
    return;
end
if ~isa(hm, 'java.util.HashMap')
    warning ('expecting a hash map or java arraylist of hashmaps');
    obj = [];
    return;
end

keys = hm.keySet;
it = keys.iterator;
for j = 1:keys.size();
    try
        fn = it.next();
        val = hm.get(fn);
        if (isjava(val))
            obj.(fn) = objFromHashmap(val);
        else
            obj.(fn) = val;
        end
    catch me
        'parse error'
        me.getReport;
    end
end
