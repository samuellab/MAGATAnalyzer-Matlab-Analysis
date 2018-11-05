function mtp2 = unCameraCalibrate (mtp, camcalinfo)
    mtp2 = mtp;
    if (isempty(camcalinfo) || ~isa(camcalinfo, 'CameraCalibration'))
        return;
    end
    
    mtp2.head = camcalinfo.camPtsFromRealPts(mtp2.head);
    mtp2.tail = camcalinfo.camPtsFromRealPts(mtp2.tail);
    mtp2.mid = camcalinfo.camPtsFromRealPts(mtp2.mid);
    mtp2.contour = camcalinfo.camPtsFromRealPts(mtp2.contour);
    mtp2.spine = camcalinfo.camPtsFromRealPts(mtp2.spine);
    mtp2.loc = camcalinfo.camPtsFromRealPts(mtp2.loc);
    
    fn = fieldnames(mtp2);
    for j = 1:length(fn)
        mtp2.(fn{j}) = double(mtp2.(fn{j}));
    end