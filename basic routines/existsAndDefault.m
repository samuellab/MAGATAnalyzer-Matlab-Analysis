function ex = existsAndDefault(varname, defaultvalue)
%function ex = existsAndDefault(varname, defaultvalue)
%
%returns true if a variable called varname exists in caller and is not
%empty
%if variable does not exist or is empty, returns false and sets the
%variable to defaultvalue (if provided) in the caller workspace
%if defaultvalue is [], this is still assigned to varname

if (evalin('caller', ['exist(''' varname ''',''var'')']) && ~isempty(evalin('caller', varname)))
    ex = true;
else
    if (exist('defaultvalue','var'))
        assignin('caller', varname, defaultvalue)
    end
    ex = false;
end
