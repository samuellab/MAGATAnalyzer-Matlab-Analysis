function [type, valid] = getJavaPtType(fname)

%     jEx = javaObjectEDT('TrackExtractionJava.Experiment', fname);
%     code = jEx.getTypeCode();
    
    code = javaMethod('getPointType','TrackExtractionJava.Experiment', fname);
    
    valid = true;
    switch (code)
        case 3
            type = LarvaTrackPoint();
        case 2
            type = MaggotTrackPoint();
        case 1
            type = ImTrackPoint();
        case 0
            type = TrackPoint();
        otherwise
            type = -1;
            valid = false;
    end
    
    
end