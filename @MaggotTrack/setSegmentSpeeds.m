function setSegmentSpeeds (track, mso)
% automatically sets run start/end speed thresholds based on speed near
% points of high curvature
% 
% function setSegmentSpeeds (track, mso)
%
% outputs: none
% inputs: 
%    track < MaggotTrack
%    mso < MaggotSegmentationOptions
%       If segmentation options (mso) are given, use those.
%       Otherwise, default to existing track.so options.  

existsAndDefault('mso', track.so);
if (mso.autoset_curv_cut)
    mso.curv_cut = mso.autoset_curv_cut_mult / median(track.getDerivedQuantity('spineLength')); %autoset_curv_cut_mult default is 5
end
track.so = mso;


% Extract all the speeds and curvatures from the track (or is it all
% tracks?)
sp = track.getDerivedQuantity(mso.speed_field);
cv = track.getDerivedQuantity('curv');

% Keep only the curvatures that are above the segmentation threshold.
highcurv = abs(cv) > mso.curv_cut;
% Issue warnings if there are zero or few points above threshold.  
if (isempty(highcurv))
    disp(['locInFile = ' num2str(track.locInFile) ' no high curvature']);
    return
end
if (sum(highcurv) < 4)
    disp(['trackNum = ' num2str(track.trackNum) '  npts = ' num2str(track.npts) '  locInFile = ' num2str(track.locInFile) ' few high curvature points']);
end

% Calculate average speed (and standard deviation) of larvae with
% curvatures above the threshold.  Then set the "stop" speed cutoff to be
% the mean + stdev.  
u = mean(sp(highcurv));
s = std(sp(highcurv));
track.so.stop_speed_cut = u + s;

% Not sure what this does.  Expands which points are above the curvature
% threshold, then cuts out the original points above threshold.  The result
% being a set of points that occur near when the larval curvature is above
% threshold?  Then take the average speed (and stdev) of these points,
% using that as the cutoff for restarting a run?  
nearhc = imdilate(highcurv, ones(5)) &~highcurv;
u = mean(sp(nearhc));
s = std(sp(nearhc));
track.so.start_speed_cut = u + s;


