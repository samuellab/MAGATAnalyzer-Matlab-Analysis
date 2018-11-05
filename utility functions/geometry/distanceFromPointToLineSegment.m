function d = distanceFromPointToLineSegment(x1, x2, p)
%function d = distanceFromPointToLineSegment(x1, x2, p)
%finds the distance from the point(s) p to the line segment defined by 
%x1, x2.  d > 0 if p is to the right of the line from x1 to x2
%
%p is a 2xN list of points

ls = sum((x2-x1).^2);
v = x2-x1;
vp = [v(2);-v(1)]; %perpendicular to the right

d = vp(1)*(p(1,:) - x1(1)) + vp(2)*(p(2,:)-x1(2)); %oriented distance to infinite line

u = (v(1)*(p(1,:) - x1(1)) + v(2)*(p(2,:)-x1(2)))/ls;
d(u < 0) = sign(d(u<0)).*sqrt((p(1,u<0)-x1(1)).^2 + (p(2,u<0)-x1(2)).^2);
d(u > 1) = sign(d(u>1)).*sqrt((p(1,u>1)-x2(1)).^2 + (p(2,u>1)-x2(2)).^2);


%{
from: http://stackoverflow.com/questions/849211/shortest-distance-between-a-point-and-a-line-segment
float minimum_distance(vec2 v, vec2 w, vec2 p) {
  // Return minimum distance between line segment vw and point p
  const float l2 = length_squared(v, w);  // i.e. |w-v|^2 -  avoid a sqrt
  if (l2 == 0.0) return distance(p, v);   // v == w case
  // Consider the line extending the segment, parameterized as v + t (w - v).
  // We find projection of point p onto the line. 
  // It falls where t = [(p-v) . (w-v)] / |w-v|^2
  const float t = dot(p - v, w - v) / l2;
  if (t < 0.0) return distance(p, v);       // Beyond the 'v' end of the segment
  else if (t > 1.0) return distance(p, w);  // Beyond the 'w' end of the segment
  const vec2 projection = v + t * (w - v);  // Projection falls on the segment
  return distance(p, projection);
}
%}