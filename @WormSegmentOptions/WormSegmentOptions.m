classdef WormSegmentOptions
    %options for segmenting worm tracks
    %tracks are segmented into periods of forward motion and periods of
    %reorientation
    %
    %first sharp turns are flagged as any point at which the velocity
    %direction is changing faster than DTHETATHRESH
    %any sharp turns that are less than JOINSTPTS apart are merged into a
    %single sharp turn
    %the velocity angle into a sharp turn is the angle PTBUFFER points before
    %the turn; and the velocity angle out is PTBUFFER points after
    %if the in and out velocity angles are within ALIGNEDTHETA, the sharp
    %turn is flagged as a reversal (if out is ~180 degrees from in) or a
    %blip (if in the same direction)
    %otherwise, flagged as an omega turn
    %all sharp turns within MINRUNTIME (in seconds) of each other are
    %grouped into a single reorientation
    %a reorientation ends at the point where dtheta/dt is less than
    %STRAIGHTTHETATHRESH after the last sharp turn in that reorientation
    properties
        dthetaHiThresh = deg2rad(60); %threshold for delta theta (rad/sec) to be considered a sharp turn
        dthetaLoThresh = deg2rad(15); %threshold for sharp turn to end
        reversalCovThresh = 1.8; %if min(covRatio) < this threshold, not a reversal, but an omega turn
        omegaCovThresh = 2.0; %if min(covRatio) > this threshold, not an omega turn -- assign to reversal or blip
        speedEndSTThresh = 0.5; %sharpturn cannot end while speed < this threshold
        joinSTpts = 3; %join sharp turns that are at most STpts - 1 apart
        ptBuffer = 3; %how many points before and after to look when computing angle into and out of turn event
        alignedTheta = deg2rad(20); %how closely aligned the before and after tracks should be to count as a reversal instead of an o-turn 
        minRunTime = 5; %minimum separation in time (in seconds) to be separate reorientations
        maxBackTime = 10; %maximum time that an animal can go backwards without reorienting -- used to correct mistakes in identifying reversals
        straightThetaThresh = deg2rad(3); %threshold for delta theta to be considered straight after reorientation (dt < thresh)    
    end
    
    methods
    end
    
end

