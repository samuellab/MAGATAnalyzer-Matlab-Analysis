classdef CheckerExperimentCalculator < handle
    %Calculates Parameters for Checker Experiments, based on provided image
    %of checkerboard and camera calibration
    %   CheckerExperimentCalculator(pictureOfProjectedCheckerboard,
    %   pictureOf1cmCheckerboard)
    %   or
    %   CheckerExperimentCalculator(pictureOfProjectedCheckerboard,
    %   viscamcalibration)
    
    properties
        cc = [];
        borderSizeInCm = 0.2; %total width of border, in cm
        rx = []; %real x axis (in cm) for morphedCheckerIm
        ry = []; %real y axis (in cm) for morphedCheckerIm
        morphedCheckerIm = []; %image of projected checkerboard, 
                               %as it would appear if camera were
                               %distortion-free and calibrated in cm
        initialCheckerIm = []; %picture of projected checkerboard, as taken by camera
        imageResolution = 0.01; %resolution of image, in cm
        globalQuantities = [];
        ai = [];
    end
    
    properties (Access = protected)
        calculated = false;
    end
    
    methods %set
        function set.borderSizeInCm(cec, value)
            if (cec.borderSizeInCm ~= value)
                cec.calculated = false; %#ok<MCSUP>
            end
            cec.borderSizeInCm = value;
        end
        function set.imageResolution(cec, value)
            if (cec.imageResolution ~= value)
                cec.calculated = false; %#ok<MCSUP>
            end
            cec.imageResolution = value;
        end
    end
    
    
    methods %constructor
        function cec = CheckerExperimentCalculator(picOfChecker, picOfCBoard)
            if isa (picOfCBoard, 'CameraCalibration')
                cec.cc = picOfCBoard;
            else
                cec.cc = CameraCalibration(picOfCBoard);
            end
            cec.initialCheckerIm = picOfChecker;
        end
        
    end
    
    methods
        calculate(cec);
        gq = getGlobals(cec, varargin);
        assignGlobals(cec, expt, varargin);
    end
    
end

