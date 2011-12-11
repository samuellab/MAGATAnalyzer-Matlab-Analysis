function [inds,why] = detectPossibleSegmentationProblems(expt)
%look for possible problems with the segmentation of various tracks
%function [inds,why] = detectPossibleSegmentationProblems(expt)
%
%we look for possible problems with the segmentation of various tracks
%signs of problems
%mean(speed) < start_speed_cut
%less than 50% of points are in runs
%mean run time is less than 50% of mean run time for experiment as a whole
%
%outputs:
%INDS: the indices of problem tracks
%WHY: a string indicating the problem found for each problem track
%inputs:
%EXPT: a member of the experiment class

k = 1;
mrt = mean(expt.gatherSubField('run', 'runTime'));
inds = [];
for j = 1:length(expt.track)
    if (mean(expt.track(j).getDerivedQuantity(expt.track(j).so.speed_field)) < expt.track(j).so.start_speed_cut)
        inds(k) = j;
        why{k} = 'slow';
        k = k+1;
        continue;
    end
    if (sum(expt.track(j).isrun ~= 0) < 0.5 * length(expt.track(j).isrun))
        inds(k) = j;
        why{k} = 'few points in run';
        k = k+1;
        continue% semicolon here?
    end
    if (mean([expt.track(j).run.runTime]) < 0.5 * mrt)
        inds(k) = j;
        why{k} = 'runs are short';
        continue;
    end
end
if (~exist('why','var'))
    why = 'you are so damn sexy, you don''t have any problems';
end

