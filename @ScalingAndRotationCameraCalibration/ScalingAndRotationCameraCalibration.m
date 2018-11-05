classdef ScalingAndRotationCameraCalibration < CameraCalibration
    %SimpleScalingCameraCalibration < handle
    %maps points from camera to real points
    %after extraction, larva position and posture are measured in pixels
    %this provides a map from each point in the image (pixel location) to 
    %physical space (usually measured in centimeters).  
    %
    %in this simplified class, you provide only a scaling factor pxpercm, which is
    %applied as cm = px/pxpercm;
    %and a rotation theta. this rotation is the rotation of the real image
    %relative to the camera image. EG -- if theta is 90, an arrow that
    %points to the right in camera space will point up in real space
    properties
        pxpercm = 1;
        theta = 0;
    end
    
    methods %constructor
        function cc = ScalingAndRotationCameraCalibration(varargin)
            %cc = CameraCalibration(checkerim)
            %cc = CameraCalibration (realx, realy, camx, camy);
            switch(length(varargin))
                case 0
                    
                case 1
                    cc.pxpercm = varargin{1};
                    
                case 2
                    cc.pxpercm = varargin{1};
                    cc.theta = deg2rad(varargin{2});
            end
            cc.P2R = [cos(cc.theta) -sin(cc.theta); sin(cc.theta) cos(cc.theta)]/cc.pxpercm;
            cc.R2P = inv(cc.P2R);
        end    
    end
    properties (SetAccess = protected)
        P2R = eye(2);
        R2P = eye(2);
    end
    methods
        function [realim,realxaxis,realyaxis] = morphCamToReal(cc, camim, varargin)
            %NOT implemented yet
            disp ('not implemented yet');
            realim = camim;
            realxaxis = (1:size(camim,2))/cc.pxpercm;
            realyaxis = (1:size(camim,1))/cc.pxpercm;
        end
        function camim = morphRealToCam(cc, realim, varargin)
            %NOT implemented yet
            disp ('not implemented yet');
            camim = realim;
        end
        
        function campts = camPtsFromRealPts(cc, realpts) 
            campts = double(cc.R2P*realpts);
        end
        function realpts = realPtsFromCamPts(cc, campts)
           realpts = double(cc.P2R*campts); 
        end
        function magfactor = pixelsPerRealUnit (cc)
            magfactor = cc.pxpercm;
        end
        function magfactor = realUnitsPerPixel (cc)
            magfactor = 1/cc.pxpercm;
        end
        
    end
end

