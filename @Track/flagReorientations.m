function flagReorientations (track, wormSegmentOptions, varargin)
% groups sharp turns into periods of reorientation
% function flagReorientations (track, wormSegmentOptions)
%
% called by SegmentTrack; not usually used directly by end user
%
% outputs: none
% inputs: 
%   TRACK: member of the track class (note MaggotTrack segmentation
%       proceeds differently)
%   WORMSEGMENTOPTIONS: segmentation options, see WormSegmentOptions class
%       for details
%  optional: 
% 'UseExistingSharpTurns', [false]/true, if true, don't find new
%   sharp turns -- so that when we load them from file, we don't nuke
%   existing sharp turns -- if false, we overwrite existing turns, unless
%   user has already started flagging
%
% 'OverwriteExistingSharpTurns', [false]/true, if true, we
%   overwrite existing turns EVEN IF USER HAS ALREADY FLAGGED THEM,
%   eliminating user codes


existsAndDefault('wormSegmentOptions', track.so);
wso = wormSegmentOptions;
track.so = wso;
UseExistingSharpTurns = false;
OverwriteExistingSharpTurns = false;
varargin = assignApplicable(varargin);


track.calculateDerivedQuantity({'eti', 'sloc', 'theta', 'deltatheta', 'ddtheta'});
if ((OverwriteExistingSharpTurns || ~any(isfinite([track.sharpTurn.userCode]))) && ~UseExistingSharpTurns)
    flagOmegaTurnsAndReversals (track, wso, varargin{:});
end
reoind = 0;
st = track.sharpTurn;
%group reorientations
%if two sharp turns are less than minRunTime apart (in seconds), they're
%part of the same reorientation
eti = track.getDerivedQuantity('eti');
inreverse = false;
for j = 1:length(st)
   switch st(j).typeCode
       case {-1,0}
           if (inreverse && eti(st(j).startInd) - eti(st(j-1).endInd) < wso.maxBackTime)
               nextreo = false;
           else
              nextreo = (j == 1 || (eti(st(j).startInd) - eti(st(j-1).endInd) > wso.minRunTime));
           end
           inreverse = false;
       case {1,2}
           if (inreverse && eti(st(j).startInd) - eti(st(j-1).endInd) < wso.maxBackTime)
               st(j).typeCode = 2;
               nextreo = false;
               inreverse = false;
           else
               nextreo = (j == 1 || (eti(st(j).startInd) - eti(st(j-1).endInd) > wso.minRunTime));
               inreverse = true;
           end
       otherwise
           j
           st(j)
           st(j).type
           pause
   end
   if (nextreo) % change from difference in central ind to difference from end to start
       if (reoind > 0)
           reo(reoind).calculateMetrics;
       end
       reoind = reoind+1;
       reo(reoind) = WormReorientation(); %#ok<AGROW>
   end
   reo(reoind).addTurn(st(j));   
   reo(reoind).calculateMetrics();
end

track.reorientation = reo;