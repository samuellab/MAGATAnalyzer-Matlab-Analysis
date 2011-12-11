function result = executeTrackFunction(eset, func, varargin)
%function result = executeTrackFunction(eset, func, varargin)
%
%result{j} = eset.expt(j).executeTrackFunction(func, varargin)
%
if (length(eset) > 1)
    if (nargout > 0 && objNargout (eset(1).expt(1).track(1), func) ~= 0)
        for j = 1:length(eset)
            [result{1:nargout}] = eset(j).executeTrackFunction(func, varargin{:});
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
        for j = 1:length(eset.expt)
            eset(j).executeTrackFunction(func, varargin{:});
        end
    end
    return;
end

if (nargout > 0 && objNargout (eset.expt(1).track(1), func) ~= 0)
    for j = 1:length(eset.expt)
        [result{1:nargout}] = eset.expt(j).executeTrackFunction(func, varargin{:});
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
    if (true || matlabpool('size') == 0)
        for j = 1:length(eset.expt)
            eset.expt(j).executeTrackFunction(func, varargin{:});
        end
    else
        parfor j = 1:length(eset.expt)
            eset.expt(j).executeTrackFunction(func, varargin{:});
        end
    end
end