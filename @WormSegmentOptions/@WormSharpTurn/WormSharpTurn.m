classdef WormSharpTurn < TrackPart
    %when a worm track makes a sharp turn, e.g. omega turn or reversal
    
    properties
         typeCode = 0; % -1 = 'omega turn'; 0 = 'double reverse or blip'; 1 = 'reversal'; 2 = 'second reversal';
         loc = [NaN;NaN];
         centralInd = 0;
         thetaIn = NaN;
         thetaOut = NaN;
         dTheta = NaN;
         userCode = NaN; %type code assigned by user
    end
    
    methods
        calculateMetrics(st);
        typeString = type(st);
        [sym,color] = symbol(st);
        userString = usertype(st);
        setUserType (st, charCode);
    end
    
    methods %constructor
        function st = WormSharpTurn(varargin)
            %st = WormSharpTurn(track, startInd, endInd)
            arglist = {'track', 'startInd', 'endInd'};
            if (nargin > 0)
                for j = 1:min(nargin, length(arglist))
                    st.(arglist{j}) = varargin{j};
                end
                if (nargin >= length(arglist))
                    %got everything we wanted, so calculate away
                    st.calculateMetrics();
                end
            end
        end
    end
    
end

