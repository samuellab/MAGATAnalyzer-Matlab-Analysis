function prettyMovie (tp, varargin) 
%function prettyMovie (tp, varargin) 
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

tp.myMovie('prettyMovie', varargin{:});
