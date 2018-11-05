classdef SimpleScalingCameraCalibration < CameraCalibration
    %SimpleScalingCameraCalibration < handle
    %maps points from camera to real points
    %after extraction, larva position and posture are measured in pixels
    %this provides a map from each point in the image (pixel location) to 
    %physical space (usually measured in centimeters).  
    %
    %in this simplified class, you provide only a scaling factor pxpercm, which is
    %applied as cm = px/pxpercm;
    properties
        pxpercm = 1;
    end
    
    methods %constructor
        function cc = SimpleScalingCameraCalibration(varargin)
            %cc = CameraCalibration(checkerim)
            %cc = CameraCalibration (realx, realy, camx, camy);
            switch(length(varargin))
                case 0
                    return;
                case 1
                    cc.pxpercm = varargin{1};
                    return;
            end
                             
        end    
    end
    
    methods
        function [realim,realxaxis,realyaxis] = morphCamToReal(cc, camim, varargin)
            realim = camim;
            realxaxis = (1:size(camim,2))/cc.pxpercm;
            realyaxis = (1:size(camim,1))/cc.pxpercm;
        end
        function camim = morphRealToCam(cc, realim, varargin)
            camim = realim;
        end
        
        function campts = camPtsFromRealPts(cc, realpts) 
            campts = realpts*cc.pxpercm;
        end
        function realpts = realPtsFromCamPts(cc, campts)
           realpts = campts/cc.pxpercm; 
        end
        function magfactor = pixelsPerRealUnit (cc)
            magfactor = cc.pxpercm;
        end
        function magfactor = realUnitsPerPixel (cc)
            magfactor = 1/cc.pxpercm;
        end
        
    end
end

