function tp = fromFile (tp, fid, loadIm, loadContour, camcalinfo)
%trackpoint.fromFile
%function tp = fromFile (fid, loadIm, loadContour, camcalinfo)
%
%loadIm, loadContour: ignored
%
%ts = tic;
try 
     if (~exist('camcalinfo', 'var'))
         camcalinfo = [];
     end

     intType = 'int32';
     floatType = 'float32';

     tp.locInFile = ftell(fid);
     tp.ind = int16(fread(fid, 1, intType));
     tp.loc = single(readPointsFromFile (fid, 1, floatType, camcalinfo));
     tp.area = single(fread(fid, 1, floatType));
     cov = (fread(fid, 4, floatType));
     tp.cov = [cov(1);cov(2);cov(4)];
catch me
    disp (me.getReport);
    disp (['fid = ' num2str(fid)]);
    disp (['current location = ' num2str(ftell(fid))]);
    disp (['point location = ' num2str(tp.locInFile)]);
    disp (['fid error message = ' ferror(fid)]);
end
    
 %{
 if all(isfinite(cov))
     [V,D] = eig(reshape(cov,[2 2]));
    if (D(1,1) > D(2,2))
        tp.cov = [D(1,1);D(2,2);atan2(V(2,1), V(1,1))];
    else
        tp.cov = [D(2,2);D(1,1);atan2(V(2,2), V(1,2))];
    end
    %{
    tp.cov = [D(2,2);D(1,1);acos(V(1,2))];
    if (tp.cov(1) < tp.cov(2))
         tp.cov = [D(1,1);D(2,2);acos(V(1,1))];
    end
    %}
 end
 %}
%toc(ts)
