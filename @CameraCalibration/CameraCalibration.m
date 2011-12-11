classdef CameraCalibration < handle
    %CameraCalibration < handle
    %maps points from camera to real points
    %after extraction, larva position and posture are measured in pixels
    %this provides a map from each point in the image (pixel location) to 
    %physical space (usually measured in centimeters).  
    %
    %either provide a list of real points and corresponding pixel locations
    %cc = CameraCalibration (realx, realy, camx, camy)
    % or provide an image of a checkerboard with squares 1 cm on a side
    %cc = CameraCalibration (im, . . .)
    % for a list of optional arguments, help calibrateCheckerboard
    
    properties
        realx;
        realy;
        camx;
        camy;
    end
    methods %constructor
        function cc = CameraCalibration(varargin)
            %cc = CameraCalibration(checkerim)
            %cc = CameraCalibration (realx, realy, camx, camy);
            switch(length(varargin))
                case 0
                    return;
                case 1
                    
                    if (isa(varargin{1}, 'CameraCalibration'))
                        cc = CameraCalibration(varargin{1}.realx, varargin{1}.realy, varargin{1}.camx, varargin{1}.camy);
                    else
                        %assume it's an image of a checkerboard
                        [realx, realy, camx, camy] = calibrateCheckerboard(varargin{1}); %#ok<*PROP>
                        cc = CameraCalibration(realx, realy, camx, camy);
                    end
                case 4
                    cc.realx = varargin{1};
                    cc.realy = varargin{2};
                    cc.camx = varargin{3};
                    cc.camy = varargin{4};
                    cc.setTSI();
                otherwise
                    try
                        [realx, realy, camx, camy] = calibrateCheckerboard(varargin{:});
                        cc = CameraCalibration(realx, realy, camx, camy);
                    catch me
                        disp(me.getReport);
                        disp('CameraCalibration(im) or CameraCalibration(realx, realy, camx, camy)');
                        cc = [];
                    end
                    return;
            end
        end    
    end
    
    methods
        [realim,realxaxis,realyaxis] = morphCamToReal(cc, camim, varargin); %take a camera picture and morph it to remove lens distortion & scale to actual size according to axes
        camim = morphRealToCam(cc, realim, realxaxis, realyaxis, varargin); %take an image defined in real space & morph it to appear as it would if the camera had taken a picture of it
        campts = camPtsFromRealPts(cc, realpts); %map real locations to a points on the camera sensor
        realpts = realPtsFromCamPts(cc, campts); %map points on the camera sensor to real locations
        magfactor = pixelsPerRealUnit (cc); %how many pixels in a cm
        magfactor = realUnitsPerPixel (cc); %how many cm/pixel
        setTSI(cc); %set tri scattered interpolation
    end
    
    properties (SetAccess = protected)
        %tri scattered interpolants
        c2rX; 
        c2rY;
        r2cX;
        r2cY;
    end
end

