function yout = interpLogicalIm(xin, xdata, ydata) 
%function yout = interpLogicalIm(xin, xdata, ydata) 


yout = logical(round(interp2(double(xdata.x), double(xdata.y), double(bwunpack(ydata)), double(xin(1,:)), double(xin(2,:)), '*linear')));
