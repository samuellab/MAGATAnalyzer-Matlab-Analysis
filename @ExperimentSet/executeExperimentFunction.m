function varargout = executeExperimentFunction(eset, func, varargin)
% executes a given method from the Experiment class on each expt
% function result = executeExperimentFunction(eset, func, varargin)
%
% func is the string name of a function in experiment class
% result{j} = eset.expt(j).(func)(varargin{:}) 
%
% outputs:
% VARARGOUT: up to the number of arguments returned by func
% in pseudo-code, [varargout{1}, varargout{2}, etc.] =
%                 [eset.expt.func(varargin)]
% inputs: 
% ESET: a member of the ExperimentSet class
% FUNC: a string that is the name of a method in the Experiment class
% VARARGIN: passed to FUNC

if (nargout > 0 && objNargout(eset.expt(1), func) ~= 0)
    
    for j = 1:length(eset.expt)   
        
         [result{1:nargout}] = eset.expt(j).(func)(varargin{:});
         allresult{j} = result;
    end
    
    for k = 1:nargout
        clear temp;
        for j = 1:length(allresult)
            temp{j} = allresult{j}{k};
        end
        varargout{k} = temp;
    end
else
    if (true || matlabpool('size') == 0)
        for j = 1:length(eset.expt)
            eset.expt(j).(func)(varargin{:});
        end
    else
        parfor j = 1:length(eset.expt)
            eset.expt(j).(func)(varargin{:});
        end
    end
end

