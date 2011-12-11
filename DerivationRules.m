classdef DerivationRules
    %Rules for Deriving Quantities in Tracks
    %no methods really
    
    %all times in seconds
    
    properties
        interpTime = 0.25; %all location data, etc. is resampled at this time scale
        smoothTime = 0.5; %sigma for gaussian smoothing filter
        derivTime = 0.25; %sigma for derivative filter (derivatives are taken of smoothed quantities)
    end
    
    methods
    end
    
end

