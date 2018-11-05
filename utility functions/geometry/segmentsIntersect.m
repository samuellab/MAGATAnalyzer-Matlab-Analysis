function si = segmentsIntersect(pt1,pt2,pt3,pt4)
%function si = segmentsIntersect(pt1,pt2,pt3,pt4)
%
%returns true if the segments defined by pt1-pt2 and pt3-pt4 intersect

d = (pt4(2,:)-pt3(2,:)).*(pt2(1,:)-pt1(1,:)) - (pt4(1,:)-pt3(1,:)).*(pt2(2,:)-pt1(2,:));
n1 = (pt4(1,:)-pt3(1,:)).*(pt1(2,:)-pt3(2,:)) - (pt4(2,:)-pt3(2,:)).*(pt1(1,:)-pt3(1,:));
n2 = (pt2(1,:)-pt1(1,:)).*(pt1(2,:)-pt3(2,:)) - (pt2(2,:)-pt1(2,:)).*(pt1(1,:)-pt3(1,:));

si = (n1./d <= 1) & (n1./d >= 0) & (n2./d <= 1) & (n2./d >= 0);
si(d == 0) = false;
inds = d == 0 & n1 == 0 & n2 == 0;
if (any(inds))
    r = pt2-pt1;
    d1(1,:) = pt3(1,:)-pt1(1,:);d1(2,:) = pt3(2,:)-pt1(2,:);
    d2(1,:) = pt4(1,:)-pt1(1,:);d2(2,:) = pt4(2,:)-pt1(2,:);
    si(inds) = (dot(r(:,inds),d1(:,inds)) <= dot(r(:,inds),r(:,inds)) & dot(r(:,inds),d1(:,inds))>0) | ...
        (dot(r(:,inds),d2(:,inds)) <= dot(r(:,inds),r(:,inds)) & dot(r(:,inds),d2(:,inds))>0);
end

% if (any(si))
%     figure(10)
%     if (size(pt1,2) > 1)
%         ind1 = find(si);
%     else
%         ind1 = 1;
%     end
%     if (size(pt3,2) > 1)
%         ind2 = find(si);
%     else
%         ind2 = 1;
%     end
%     plot ([pt1(1,ind1) pt2(1,ind1)], [pt1(2,ind1) pt2(2,ind1)], 'b-', [pt3(1,ind2) pt4(1,ind2)], [pt3(2,ind2) pt4(2,ind2)], 'r-', ...
%         pt1(1,ind1), pt1(2,ind1), 'bo', pt3(1,ind2), pt3(2,ind2), 'ro');
%     d(si)
%     n1(si)
%     n2(si)
% end