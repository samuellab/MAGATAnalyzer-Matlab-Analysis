function draw(headSwing, varargin)

if (headSwing.accepted)
    c = 'g-';
else
    c = 'r-';
end
plot (headSwing.track.dq.shead(1,headSwing.inds), headSwing.track.dq.shead(2,headSwing.inds), c);