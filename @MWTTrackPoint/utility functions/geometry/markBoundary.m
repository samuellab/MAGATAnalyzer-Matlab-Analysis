function im = markBoundary (im, x0, y0, x1, y1)
if (length(x0) > 1)
    if (nargin == 3)
        for j = 1:(length(x0)-1)
            im = markBoundary (im, x0(j), y0(j), x0(j+1), y0(j+1));
        end
    else
        for j = 1:(length(x0))
            im = markBoundary (im, x0(j), y0(j), x1(j), y1(j));
        end
    end
    return;
end

dx = abs(x1 - x0);
dy = abs(y1 - y0);
x = floor(x0);
y = floor(y0);
n = 1;
if (dx == 0)
    x_inc = 0;
    err = Inf;
else
    if (x1 > x0)
        x_inc = 1;
        n = n + floor(x1) - x;
        err = (floor(x0) + 1 - x0)*dy;
    else
        x_inc = -1;
        n = n + x - floor(x1);
        err = (x0 - floor(x0)) * dy;
    end
end

if (dy == 0)
    y_inc = 0;
    err = - Inf;
else
    if (y1 > y0)
        y_inc = 1;
        n = n + floor(y1) - y;
        err = err - (floor(y0) + 1 - y0) * dx;
    else
        y_inc = -1;
        n = n + y - floor(y1);
        err = err - (y0 - floor(y0)) * dx;
    end
end
for m = 1:n
  %  figure(1);
  %  pcolor(im); axis equal; hold on
      if (y >= 1 && y <= size(im,1) && x >= 1 && x <= size(im,2))
        im(y,x) = divideSquareByLine([x;y], [x+1;y+1], [x0;y0], [x1;y1]); 
      end
 %   title ([num2str(x) ',' num2str(y) '; ' num2str(x0) ',' num2str(y0) ' - ' num2str(x1) ',' num2str(y1)]);
 %   pause;
    if (err > 0)    
        y = y + y_inc;
        err = err - dx;
    else
        x = x + x_inc;
        err = err + dy;
    end
end


%{
adaptation of:
#include <limits> // for infinity
void raytrace(double x0, double y0, double x1, double y1)
{
    double dx = fabs(x1 - x0);
    double dy = fabs(y1 - y0);

    int x = int(floor(x0));
    int y = int(floor(y0));

    int n = 1;
    int x_inc, y_inc;
    double error;

    if (dx == 0)
    {
        x_inc = 0;
        error = std::numeric_limits<double>::infinity();
    }
    else if (x1 > x0)
    {
        x_inc = 1;
        n += int(floor(x1)) - x;
        error = (floor(x0) + 1 - x0) * dy;
    }
    else
    {
        x_inc = -1;
        n += x - int(floor(x1));
        error = (x0 - floor(x0)) * dy;
    }

    if (dy == 0)
    {
        y_inc = 0;
        error -= std::numeric_limits<double>::infinity();
    }
    else if (y1 > y0)
    {
        y_inc = 1;
        n += int(floor(y1)) - y;
        error -= (floor(y0) + 1 - y0) * dx;
    }
    else
    {
        y_inc = -1;
        n += y - int(floor(y1));
        error -= (y0 - floor(y0)) * dx;
    }

    for (; n > 0; --n)
    {
        visit(x, y);

        if (error > 0)
        {
            y += y_inc;
            error -= dx;
        }
        else
        {
            x += x_inc;
            error += dy;
        }
    }
}
from http://playtechs.blogspot.com/2007/03/raytracing-on-grid.html
%}

