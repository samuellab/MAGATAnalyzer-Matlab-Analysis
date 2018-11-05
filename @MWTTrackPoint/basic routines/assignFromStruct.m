function y = assignFromStruct(x)
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
fn = fieldnames(x);

for i = 1:length(fn)
    if (any(strcmp(fn{i},vars)))
        assignin('caller', fn{i}, x.(fn{i}));
        used = [used i]; %#ok<AGROW>
    end
    
end
if (~isempty(used))
    inds = setdiff(1:length(fn), used);
    for j = inds
        y.(fn{inds{j}}) = x.(fn{inds(j)});
    end
end

%############################################

