function [a, b] = minbrak(fn, a, b)
%function [a, b] = minbrak(fn, a, b)
%
%brackets the minimum of fn; adapted from numerical recipes in c
GOLD = 1.62;
fa = fn(a);
fb = fn(b);

if (fa < fb)
    c = a;
    a = b;
    b = c;
    fc = fa;
    fa = fb;
    fb = fc;
end

c = (1+GOLD)*b - GOLD*a;
fc = fn(c);
while (fc < fb) 
    r = (b-a)*(fb-fc);
    q = (b-c)*(fb-fa);
    u = b - ((b-c)*q - (b-a)*r)/2*(q-r + eps*sign(q-r));
    ulim = b + 10*c-b;
    if ((b-u)*(u-c) > 0)
        fu = fn(u);
        if (fu < fc)
            a = b;
            b = c;
            return;
        end;
        if (fu > fb)
            b = u;
            return;
        end
        u = c + GOLD*(c-b);
        fu = fn(u);
    else
        if ((c - u) *(u - ulim) > 0)
            fu = fn(u);
            if (fu < fc)
                b = c;
                c = u;
                u = GOLD*(c-b);
                fb = fc;
                fc = fu;
                fu = fn(u);
            end
        else
            if ((u - ulim)*(ulim - c) > 0)
                u = ulim;
                fu = fn(u);
            else
                u = c + GOLD*(c-b);
                fu = fn(u);
            end
        end
    end
    a = b;
    b = c;
    c = u;
    fa = fb;
    fb = fc;
    fc = fu;
end
b = c;
