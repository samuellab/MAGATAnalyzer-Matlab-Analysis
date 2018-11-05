function x = assignApplicableAndFixTypes(x)
% PVPMOD             - evaluate parameter/value pairs
% pvpmod(x) assigns the value x(i+1) to the parameter defined by the
% string x(i) in the calling workspace if and only if the calling function
% already has a parameter by the name x{i} defined
% otherwise, does nothing
% unused parameter-value pairs are returned in x
% This is useful to evaluate 
% <varargin> contents in an mfile, e.g. to change default settings 
% of any variable initialized before pvpmod(x) is called.
%
% modified by marc gershow from pvpmod by
% (c) U. Egert 1998

%############################################
% this loop is assigns the parameter/value pairs in x to the calling
% workspace.
used = [];
vars =  evalin('caller', 'who');
%s = evalin('caller', 'whos');
%vars = {s.name}';
if ~isempty(x)
    skipnext = false;
   for i = 1:size(x,2)
       if skipnext
           skipnext = false;
           continue;
       end
       if (ischar(x{i}) && any(strcmp(x{i},vars)))   
 %         ind = find(strcmp(x{i}, vars), 1, 'first');
          oldval = evalin('caller', x{i});
          val = fixType(x{i+1}, oldval);
          assignin('caller', x{i}, val);
          used = [used i];
          skipnext = true;
       end
   end;
end;
if (~isempty(used))
    used = [used used+1];
    inds = setdiff(1:length(x), used);
    x = x(inds);
end

%############################################

function val = fixType(val, oldval)
if ((isnumeric(val) && isnumeric(oldval)) || isa(val, class(oldval)))
    return;
end
if(ischar(val))  
    if (isnumeric(oldval))
        val = str2double(val);
        return;
    end
    if (islogical(oldval))
        if (strcmpi(val, 'true'))
            val = true;
            return;
        end
        if (strcmpi(val, 'false'))
            val = false;
            return;
        end
        try
            v = val;
            val = logical(str2double(val));
        catch
            disp (['could not update logical value to ' v]);
            val = oldval;
            return;
        end
        return;
    end
    return
end

if (islogical(oldval))
    v = val;
    try
        val = logical(val);
    catch
        disp (['could not update logical value to ' v]);
        val = oldval;
        return
    end
end
            
