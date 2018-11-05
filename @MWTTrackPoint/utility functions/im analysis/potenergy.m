function upot = potenergy(ctr, energyim)
%function upot = potenergy(ctr, energyim)

upot = sum(interp2(energyim, ctr(1,:), ctr(2,:), '*linear', 100*max(abs(energyim(:)))));
