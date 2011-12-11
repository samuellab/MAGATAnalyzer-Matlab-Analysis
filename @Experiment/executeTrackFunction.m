function varargout = executeTrackFunction(expt, func, varargin)
%loops through all tracks and executes func(expt.track(j), varargin)
%function varargout = executeTrackFunction(expt, func, varargin)
%
%inputs:
%EXPT: a member of the Experiment class
%FUNC: the function to execute
%   if func is a string, calls
%   expt.track(j).(func)(varargin{:})
%   if func is not a string, it is assumed to be a function handle; call is
%   func(expt.track(j),varargin{:})
%VARARGIN: additional arguments passed to func
%outputs: varargout
%   varargout is a cell array of cells. 
%   varargout{1} is a compilation of the first return output, etc.
%   in pseduo-code 
%    [varargout{1}{j}, varargout{2}{j},...] = expt.track(j).(func)(varargin{:})

if ischar(func)
    if (nargout > 0 && ~isempty(expt.track) && objNargout (expt.track(1), func) ~= 0)
        for j = 1:length(expt.track)
            [result{1:nargout}] = expt.track(j).(func)(varargin{:});
             allresult{j} = result;
        end
        for k = 1:nargout
            clear temp;
            for j = 1:length(allresult)
                temp(j) = allresult{j}{k};
            end
            varargout{k} = temp;
        end
    else
        for j = 1:length(expt.track)
            expt.track(j).(func)(varargin{:});
        end
    end
else
    if (nargout > 0 && nargout(func) ~= 0)
        for j = 1:length(expt.track)
            [result{1:nargout}]  = func(expt.track(j),varargin{:});
             allresult{j} = result;
        end
        for k = 1:nargout
            clear temp;
            for j = 1:length(allresult)
                temp(j) = allresult{j}{k};
            end
            varargout{k} = temp;
        end
    else
        for j = 1:length(expt.track)
           func(expt.track(j),varargin{:});
        end
    end
end
        