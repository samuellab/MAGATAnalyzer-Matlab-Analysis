function varargout = nameInList (list, varargin)
% varargout = Track.nameInList (list, varargin)
% helper function for validDQName

switch (nargin-1)
    case 0
        varargout{1} = list;
    case 1
        if (iscell(varargin{1}))
            tfarray = repmat(false, size(varargin{1}));
            for j = 1:length(varargin{1})
                tfarray(j) = Track.nameInList(list, varargin{1}{j});
            end
            varargout{1} = tfarray;
        else 
            varargout{1} = (ischar(varargin{1}) && (any(strcmp(varargin{1}, list))));
        end
    otherwise
        for j = 1:(nargin-1)
            varargout{j} = Track.nameInList(list, varargin{j});
        end
end