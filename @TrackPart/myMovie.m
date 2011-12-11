function myMovie (tp, mtype, varargin) 
%function myMovie (tp, mtype, varargin) 
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

nopause = false;
showind = false;
varargin = assignApplicable(varargin);

if (length(tp) > 1)
    for j = 1:length(tp)
        if (showind)
            disp(num2str(j));
        end
        tp(j).myMovie(mtype, varargin{:});
        if (~nopause)
            pause;
        end
    end
    return;
end

nreps = 1;
ptbuffer = 10;
delayTime = 0.05;
%class(tp)
if (isa(tp, 'HeadSwing'))
    delayTime = 0.15;
    ptbuffer = 3;
end
if (isa(tp, 'MaggotReorientation'))
    delayTime = 0.1;
    ptbuffer = 3;
end

varargin = assignApplicable(varargin);
p2i = tp.track.getDerivedQuantity('mapInterpedToPts');

si = max(p2i(tp.startInd) - ptbuffer, 1);
ei = min(p2i(tp.endInd) + ptbuffer, tp.track.npts);
    
figure(gcf);
for j = 1:nreps
    tp.track.(mtype)('inds', si:ei, 'delayTime', delayTime, 'highlightInds', tp.startInd:tp.endInd, varargin{:});
end