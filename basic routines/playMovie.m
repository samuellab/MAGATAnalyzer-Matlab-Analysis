function playMovie (fn,start,stop)

blocksize = 1000;
inf = aviinfo(fn);

if (~exist('start','var') || isempty(start))
    start = 1;
end
if (~exist('stop','var') || isempty(stop))
    stop = inf.NumFrames;
end

ind1 = start;
ind2 = ind1+blocksize-1;
mov = aviread(fn,ind1:ind2);

for j = start:stop
    if (mod((j),blocksize) == 0)
        ind1 = j;
        ind2 = ind1+blocksize;
        if (ind2 > stop)
            ind2 = stop;
        end
        mov = aviread(fn,ind1:ind2);
    end
    imagesc(double(mov(mod(j-1,blocksize)+1).cdata));
    title (['Frame ' num2str(j)]);
    colormap gray
    pause (0.02);
end