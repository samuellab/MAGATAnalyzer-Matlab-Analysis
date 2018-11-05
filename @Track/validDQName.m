function varargout = validDQName (varargin)
% tells whether a name is a valid derived quantity/lists all valid DQs
% function varargout = validDQName (varargin)
%
% tf = validDQName(name) -- true or false if name is a valid derived
%       field for a track
% tflist = validDQName(namelist) -- array of true or falses, namelist is a cell
% tflist = validDQName('name1', 'name2') -- nargout = nargin t/f
% namelist = validDQName() - returns a cell containing all valid field names

validfieldlist = {'eti', 'iloc','sloc','vel', 'speed', 'nspeed', 'nvel', 'vnorm', 'speed_diff_local', 'theta','deltatheta', 'ddtheta', 'acc', 'curv',...
                    'pathLength','displacement','icov', 'covRatio', 'covTheta','scov', 'scovRatio', 'scovTheta', 'covMajor', 'adjspeed' ...
                    'covMinor', 'scovMajor', 'scovMinor','iarea', 'sarea', 'totalTime', 'lrdtheta', 'dcovRatio','xloc','yloc', 'etiFromTrackStart','mapinterpedtopts',...
                    'isrun', 'iscollision','iiscollision'};
                
                
varargout{:} = Track.nameInList(validfieldlist, varargin{:});
%{
switch nargin
    case 0
        varargout{1} = validfieldlist;
    case 1
        if (iscell(varargin{1}))
            tfarray = repmat(false, size(varargin{1}));
            for j = 1:length(varargin{1})
                tfarray(j) = Track.validDQName(varargin{1}{j});
            end
            varargout{1} = tfarray;
        else 
            varargout{1} = (ischar(varargin{1}) && (any(strcmp(varargin{1}, validfieldlist))));
        end
    otherwise
        for j = 1:nargin
            varargout{j} = Track.validDQName(varargin{j});
        end
end
%}
            
        
