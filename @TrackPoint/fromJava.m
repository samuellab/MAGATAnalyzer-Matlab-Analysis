function tp = fromJava(tp, jTP, loadIm, loadContour, camcalinfo)

    if (~exist('camcalinfo', 'var'))
        camcalinfo = [];   
    end
    
    %NOT SETTING LOC IN FILE
    tp.ind = jTP.getFrameNum;
    tp.loc = properCoords([jTP.getX; jTP.getY], camcalinfo);
    tp.area = jTP.getArea; %In pixels^2
    %NOT SETTING COV
    
    
end