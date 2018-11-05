function setTSI(cc)
%function setTSI(cc)
%
% used internally by CameraCalibration; do not call directly unless you know
% what you are doing
cc.camx = makecolumn(cc.camx);
cc.camy = makecolumn(cc.camy);
cc.realx = makecolumn(cc.realx);
cc.realy = makecolumn(cc.realy);

[cx, cy, rx, ry] = guessOutsideHull (cc.camx, cc.camy, cc.realx, cc.realy, [min(cc.realx) max(cc.realx)], [min(cc.realy) max(cc.realy)]);
cc.r2cX = TriScatteredInterp (rx, ry, cx);
cc.r2cY = TriScatteredInterp (rx, ry, cy);

%1/16 - added 10 pixels to border region to take care of weird stuff from extraction, but
%I need to check what happened during extraction
[rx, ry, cx, cy] = guessOutsideHull (cc.realx, cc.realy, cc.camx, cc.camy, [min(cc.camx)-10 max(cc.camx)+10], [min(cc.camy)-10 max(cc.camy)+10]);
cc.c2rX = TriScatteredInterp (cx, cy, rx);
cc.c2rY = TriScatteredInterp (cx, cy, ry);
%cc.c2rX = TriScatteredInterp (cc.camx, cc.camy, cc.realx);
%cc.c2rY = TriScatteredInterp (cc.camx, cc.camy, cc.realy);

