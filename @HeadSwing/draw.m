function draw(headSwing, varargin)
Colors.accepted = 'g-';
Colors.rejected = 'r-';
varargin = assignApplicable(varargin);

if (headSwing.accepted)
    c = Colors.accepted;
else
    c = Colors.rejected;
end

plot (headSwing.track.dq.shead(1,headSwing.inds), headSwing.track.dq.shead(2,headSwing.inds), c);