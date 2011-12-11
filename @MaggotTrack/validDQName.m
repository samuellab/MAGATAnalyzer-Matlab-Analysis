function varargout = validDQName (varargin)
%@Track/validDQName
%static function that tells whether a name is a valid derived quantity
%function varargout = validDQName (varargin)
%
%tf = validDQName(name) - true or false if name is a valid derived
%field for a track
%tflist = validDQName(namelist) - array of true or falses, namelist is a cell
%tflist = validDQName('name1', 'name2') - nargout = nargin t/f
%namelist = validDQName() - returns a cell containing all valid field names

persistent validfieldlist;
if (isempty(validfieldlist))
    validfieldlist = {'ihead', 'shead', 'itail', 'stail', 'imid', 'smid', 'ihtValid', 'ibodytheta', 'sbodytheta', ...
                    'vtail', 'vhead', 'vmid', 'sptail', 'sphead', 'spmid', 'vel_dp','imhdir','itmdir','smhdir','stmdir','vheadperp','spheadperp'...
                    'ispine','spineLength','spineWidth','periAmp', 'periFreq', 'periTau', ...%'periPhase', 'periMean', ...
                    'spineDist','spineCurv', 'spineTheta', 'sspineTheta', 'dspineCurv', 'dspineTheta','dbodytheta', 'dsbodytheta'};
    validfieldlist = union(validfieldlist, Track.validDQName);
end

varargout{:} = Track.nameInList(validfieldlist, varargin{:});

            
        
