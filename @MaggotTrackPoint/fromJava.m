function tp = fromJava(tp, jTP, loadIm, loadContour, camcalinfo)

    %tp = fromFile@ImTrackPoint(tp, fid, loadImage, loadContour, camcalinfo);
    tp = fromJava@ImTrackPoint(tp, jTP, loadIm, loadContour, camcalinfo);
    
    %NOT SETTING TARGETAREA
    tp.threshold = jTP.getThresh();
    tp.htValid = jTP.getHTValid();
    if (loadContour)
%         tp.contour = properCoords(JavaConPts2Array(jTP.getContour())', camcalinfo);
        tp.contour = properCoords(jTP.getContourArray(), camcalinfo);
    end
    
    tp.head = properCoords(jTP.getHead(), camcalinfo);
    tp.mid = properCoords(jTP.getMid(), camcalinfo);
    tp.tail = properCoords(jTP.getTail(), camcalinfo);
    if (tp.htValid)
        tp.spine = properCoords(jTP.getMidlineArray(), camcalinfo);
    end
%     if (tp.htValid)
%         tp.head = properCoords([jTP.getHead().getX(); jTP.getHead().getY()], camcalinfo);
%         tp.mid = properCoords([jTP.getMid.getX(); jTP.getMid.getY()], camcalinfo);
%         tp.tail = properCoords([jTP.getTail.getX(); jTP.getTail.getY()], camcalinfo);
        
        %CHANGE TO GET MIDLINE ARRAY 
%     else
%         tp.head = [NaN;NaN];
%         tp.mid = [NaN;NaN];
%         tp.tail = [NaN;NaN];
%         tp.spine = [NaN;NaN];
%     end
    
    
    
end