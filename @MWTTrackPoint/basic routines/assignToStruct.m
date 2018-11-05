function [s,x] = assignToStruct(s, x)
%function [s,x] = assignToStruct(s, x)
%
%if x{i} is a field in s, 
%s.(x{i}) = x{i+1}
%
%removes any parameter, value pairs from x and returns
% modified by marc gershow from pvpmod by
% (c) U. Egert 1998

%############################################
% this loop is assigns the parameter/value pairs in x to the calling
% workspace.
used = [];
if ~isempty(x)
    skipnext = false;
   for i = 1:size(x,2)
       if skipnext
           skipnext = false;
           continue;
       end
       if (ischar(x{i}) && isfield(s,x{i}))         
          s.(x{i}) = x{i+1};
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

