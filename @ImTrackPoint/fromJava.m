function tp = fromJava(tp, jTP, loadIm, loadContour, camcalinfo)

    %tp = fromJava@TrackPoint(tp, jTP, loadIm, loadContour, camcalinfo);
    tp = fromJava@TrackPoint(tp, jTP, loadIm, loadContour, camcalinfo);
    
    r = jTP.getRect();
    tp.imOffset = [r(1) r(2)];
    if (loadIm)
        tp.imData = jTP. getRawIm().getIntArray()';%Probably need to convert(/transpose?) this
    end
    
end