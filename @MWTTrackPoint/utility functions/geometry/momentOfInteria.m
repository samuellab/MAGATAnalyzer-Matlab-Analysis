function moi = momentOfInteria(midline, c1, c2)
%function moi = momentOfInteria(midline, c1, c2)

dt = DelaunayTri(midline');

%dt1 = DelaunayTri(c1');
%dt2 = DelaunayTri(c2');

%nn1 = dt1.nearestNeighbor(midline');
%nn2 = dt2.nearestNeighbor(midline');
nn1 = dt.nearestNeighbor(c1');
nn2 = dt.nearestNeighbor(c2');

%compute areas
v1 = c1 - midline(:,nn1);
v1 = interp1(v1', linspace(1,length(v1), length(v1)-1))';
w1 = diff(c1,[],2);
a1 = abs(v1(1,:).*w1(2,:) - v1(2,:).*w1(1,:));
d1 = sqrt(sum(v1.^2));

v2 = c2 - midline(:,nn2);
v2 = interp1(v2', linspace(1,length(v2), length(v2)-1))';
w2 = diff(c2,[],2);
a2 = abs(v2(1,:).*w2(2,:) - v2(2,:).*w2(1,:));
d2 = sqrt(sum(v2.^2));



moi = sum(a2.*d2 + a1.*d1);
