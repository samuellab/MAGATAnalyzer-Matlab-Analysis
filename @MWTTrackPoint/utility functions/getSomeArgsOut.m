function varargout = getSomeArgsOut (args, fun, varargin) %#ok<INUSD,STOUT>
%function varargout = getSomeArgsOut (args, fun, varargin)
%
%calls fun(varargin{:}) and returns the arguments selected by args
%
%example usage:
%>>getSomeArgsOut(1, @max, (1:7).^2)
%ans =
%     49
%>>getSomeArgsOut(2, @max, (1:7).^2)
%ans =
%     7

arglist = '[';
for j = 1:max(args)
    if (j ~= 1)
        arglist = [arglist ',']; %#ok<*AGROW>
    end
    if (any(args == j))
        arglist = [arglist 'arg' num2str(j)];
    else
        arglist = [arglist '~'];
    end
end
arglist = [arglist ']'];

eval([arglist ' = fun(varargin{:});']);

for j = 1:length(args)
    eval( ['varargout{j} = arg' num2str(args(j)) ';']);
end

    
    
    