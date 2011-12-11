function fixHTOrientation(track, varargin)
% legacy code; fixes extraction mistakes that swap head & tail
% @MaggotTrack/fixHTOrientation(track, varargin)
% maggotTrack.fixHTOrientation(varargin)
%
% fixes head/tail orientations;
% the following parameter/value pairs can be passed in to modify shown
% defaults
%
% outputs: none
% inputs: TRACK < MaggotTrack
% optional arguments (pass as parameter/value pairs)
%   mintime = 10 (seconds running backwards that signals a problem)
%   speedthresh = 33rd percentile of speed -- if it's going slower than this, 
%               we don't worry that it's going backwards
%   dpthresh = -0.5 if the norm of the mid-head vector dot the norm of the velocity
%         is below this threshold, it's going the wrong way
if (~isempty(track.dq))
    existingFields = fieldnames(track.dq);
else
    existingFields = {};
end

mintime = 10; %seconds minimum time for something to be bad before it needs to be fixed
pt = [track.pt];
if (sum([pt.htValid]) < 4 || sum([pt.htValid]) < 4*mintime/track.dr.interpTime)
    return; %track is too short
end
track.calculateDerivedQuantity({'eti', 'vel_dp', 'speed'});

sp = track.getDerivedQuantity('speed');
sp2 = sort(sp);
speedthresh = sp2(round(length(sp2)/3));
dpthresh = -0.5;

assignApplicable(varargin);
minpts = mintime/track.dr.interpTime;
dp = track.getDerivedQuantity('vel_dp');

inds = find(sp > speedthresh & isfinite(dp));

if (length(inds) < 2)
    disp ('no points above speed cut in fixHTOrientation');
    return;
end
dp = double(interp1(inds,dp(inds), 0:(length(dp)+1), 'nearest', 'extrap'));
dpfilt = medfilt1(dp, minpts);
badinds = dpfilt < dpthresh;
badinds = find(imdilate(badinds, ones([1 3])));
badinds = badinds((badinds >= 1) & (badinds <= length(track.dq.eti))); 
badtimes = track.dq.eti((badinds)); %time when it's going backwards

%find the original points that are in that time range
%note we are assuming that the interpolation time <= sampling time, which
%should be valid; we space out by thirds to make sure we cover all points
%note unique sorts, so there is no need to do that here
badtimes = unique([badtimes, badtimes-track.dr.interpTime/3, badtimes+track.dr.interpTime/3]);
badinds = unique(interp1([pt.et], 1:track.npts, badtimes,'nearest', 'extrap'));

%swap head and tail in bad time range
if (~isempty(badinds))
    allpts = track.pt;
    for j = badinds
        temp = allpts(j).head;
        allpts(j).head = allpts(j).tail;
        allpts(j).tail = temp;
        if (isfield(allpts(j), 'spine') && ~isempty(allpts(j).spine))
            allpts(j).spine = allpts(j).spine(end:-1:1,:);
        end
    end
    track.pt = allpts;
end

%only recalculate maggot track params (track params don't involve head or
%tail, so won't change)
fnames = setdiff(MaggotTrack.validDQName, Track.validDQName);
track.recalculateDerivedQuantities(fnames{:});

%{
%clear all fields from dq, then rederive any that existed before
track.dq = [];
if ~isempty(existingFields) 
    track.calculateDerivedQuantity(existingFields);
end
%}