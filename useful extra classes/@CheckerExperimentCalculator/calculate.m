function calculate(cec)
%function calculate(cec)

[im,rxa,rya] = cec.cc.morphCamToReal(cec.initialCheckerIm,'realresolution',cec.imageResolution);

cec.morphedCheckerIm = im;
cec.rx = rxa;
cec.ry = rya;

borderpts = cec.borderSizeInCm/cec.imageResolution;
cornerpts = borderpts*1.5;
cec.ai = analyzeCheckerIm(cec.morphedCheckerIm, borderpts, cornerpts);
cec.ai.distToBorder = cec.imageResolution * cec.ai.distToBorder;
cec.ai.distToDark = cec.imageResolution * cec.ai.distToDark;
cec.ai.distToLight = cec.imageResolution * cec.ai.distToLight;

cec.globalQuantities = cec.getGlobals;


cec.calculated = true;