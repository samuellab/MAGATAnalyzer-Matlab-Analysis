classdef MaggotSegmentOptions
    %Options for segmenting a maggot track
    
    properties
        
        curv_cut = 0.4; %if track curvature > curv_cut, end a run
        theta_cut = pi/2; %if body theta > theta_cut, end a run
        speed_field = 'speed'; %name of field that contains speed to use (normally 'speed')
        stop_speed_cut = 2; %if speed < stop_speed_cut, end a run
        start_speed_cut = 3; %if speed > start_speed_cut && vel_dp > aligned_dp, run can start
        aligned_dp = cosd(45); %if vel_dp > aligned_dp && speed > start_speed_cut, run can start
        
        minRunTime = 2.5; %if run is less than this many seconds, discard 
        minRunLength = 0; %if run is less than this distance, discard
        
    %    minStopTime = 2; %not used right now
      %  headswing_dtheta_dt = deg2rad(10); %degrees per second
        headswing_start = deg2rad(20); %if body theta > headswing_start and not in a run, start a headswing
                                       
        headswing_stop = deg2rad(10); %if body theta < headswing_stop (or changes sign), end headswing
     %   avgAngleTime = 1; %not used right now
    end
    
    methods
    end
    
end

