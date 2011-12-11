function varargout = validDQName (varargin)
%@MartaTrack/validDQName
%static function that tells whether a name is a valid derived quantity
%function varargout = validDQName (varargin)
%
%tf = validDQName(name) - true or false if name is a valid derived
%field for a track
%tflist = validDQName(namelist) - array of true or falses, namelist is a cell
%tflist = validDQName('name1', 'name2') - nargout = nargin t/f
%namelist = validDQName() - returns a cell containing all valid field names

validfieldlist = {'smoothVel', 'smoothSpeed','spmax', 'spmin', 'spamp'};
validfieldlist = union(validfieldlist, MaggotTrack.validDQName);

varargout{:} = Track.nameInList(validfieldlist, varargin{:});

            
        
