function playMovie (tp, varargin) 
%function playMovie (tp, varargin) 
%
%plays movie, or sequence of movies, for the segments of trackpart
%tp can be a single trackpart or an array of trackparts
%
%in addition to the argument pairs that can be passed to track.playMovie 
%you can additionally pass
%
%'nreps', n (default 1): number of times to repeat playing each trackpart
%movie
%'ptbuffer', n (default 10): number of points before and after inds to play
%
%'nopause', t/f (default false): if true, don't pause between trackparts if
%tp is an array

tp.myMovie('playMovie', varargin{:});

%{
nopause = false;
varargin = assignApplicable(varargin);

if (length(tp) > 1)
    for j = 1:length(tp)
        tp(j).playMovie(varargin{:});
        if (~nopause)
            pause;
        end
    end
end

nreps = 1;
ptbuffer = 10;
varargin = assignApplicable(varargin);

p2i = tp.track.getDerivedQuantity('mapInterpedToPts');

si = max(p2i(tp.startInd) - ptbuffer, 1);
ei = min(p2i(tp.endInd) + ptbuffer, tp.track.npts);
    
figure(gcf);
for j = 1:nreps
    tp.track.playMovie('inds', si:ei, varargin{:});
end
%}