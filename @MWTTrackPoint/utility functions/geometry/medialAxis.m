function pts = medialAxis(cpts)
constraints = [1:length(cpts);[2:length(cpts) 1]];
dt = DelaunayTri(double(cpts'), constraints');
inside = dt.inOutStatus();
tr = TriRep(dt(inside, :), dt.X);

% Construct a set of edges that join the circumcenters of neighboring
% triangles; the additional logic constructs a unique set of such edges.
numt = size(tr,1);
T = (1:numt)';
neigh = tr.neighbors();
cc = tr.circumcenters();
xcc = cc(:,1);
ycc = cc(:,2);
idx1 = T < neigh(:,1);
idx2 = T < neigh(:,2);
idx3 = T < neigh(:,3);
neigh = [T(idx1) neigh(idx1,1); T(idx2) neigh(idx2,2); T(idx3) neigh(idx3,3)]';

pts = cc(neigh,:)';

clf;
triplot(tr, 'g');
hold on;
plot(xcc(neigh), ycc(neigh), '-r', 'LineWidth', 1.5);