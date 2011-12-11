function segmentTracks(expt, segmentOptions)
%segments tracks into runs and reorientation by calling Track.segmentTrack
%function segmentTracks(expt, segmentOptions)
%
%ouputs: none
%inputs: 
%EXPT: a member of the Experiment class
%SEGMENTOPTIONS: (optional), the segment options to be applied to every
%   track.  If empty, the track's segmentation options (not the
%   experiment's) are used.  Thus, to use the experiment's default
%   segmentation options, expt.segmentTracks(expt.so)

if (existsAndDefault('segmentOptions', []))
    expt.so = segmentOptions;
end

expt.executeTrackFunction('segmentTrack', segmentOptions);