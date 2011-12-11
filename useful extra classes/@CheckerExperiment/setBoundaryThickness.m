function setBoundaryThickness(expt, thickness)
%function setBoundaryThickness(expt, thickness)
%sets whether or not a point is on the boundary based on the given
%thickness; boundaryDistance must have already been calculated 

expt.boundaryThickness = thickness;

corners = false(expt.imsize);
inds = sub2ind (size(corners), round(expt.structvar(:,2)), round(expt.structvar(:,1)));
corners(inds) = true;
corners = imdilate(corners, ones(round(expt.boundaryThickness/2)));

[x,y] = meshgrid(1:expt.imsize(2),1:expt.imsize(1));
xdata(:,:,1) = x;
xdata(:,:,2) = y;
for j = 1:length(expt.track)
    if (isfield(expt.track(j).dq, 'boundaryDistance'))
        %expt.track(j).dq.boundaryDistance = abs(expt.track(j).dq.boundaryDistance);
        iscorner = GlobalQuantity.twoDinterpolation(expt.track(j).getDerivedQuantity('sloc'), xdata, corners);
        expt.track(j).dq.onboundary = abs(expt.track(j).dq.boundaryDistance) <= thickness/2 & ~iscorner;
       % expt.track(j).dq.boundaryDistance(iscorner) = -expt.track(j).dq.boundaryDistance(iscorner);
    end
end