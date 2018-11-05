%LUMINOS	Photopic luminosity function
%
% 	P = LUMINOS(lam)
%
% 	Return photopic luminosity (0..1) for wavelength specified.  
%	Lam maybe a vector.
%
% lambda is in nm
%
%	Copyright (c) Peter Corke, 1999  Machine Vision Toolbox for Matlab
%   modified by marc gershow 3/2/2012
  function lu = luminos(lam)
	  tab = [
	  3.8000000e-007  0.0000000e+000
	  3.8500000e-007  1.0000000e-004
	  3.9000000e-007  1.0000000e-004
	  3.9500000e-007  2.0000000e-004
	  4.0000000e-007  4.0000000e-004
	  4.0500000e-007  6.0000000e-004
	  4.1000000e-007  1.2000000e-003
	  4.1500000e-007  2.2000000e-003
	  4.2000000e-007  4.0000000e-003
	  4.2500000e-007  7.3000000e-003
	  4.3000000e-007  1.1600000e-002
	  4.3500000e-007  1.6800000e-002
	  4.4000000e-007  2.3000000e-002
	  4.4500000e-007  2.9800000e-002
	  4.5000000e-007  3.8000000e-002
	  4.5500000e-007  4.8000000e-002
	  4.6000000e-007  6.0000000e-002
	  4.6500000e-007  7.3900000e-002
	  4.7000000e-007  9.1000000e-002
	  4.7500000e-007  1.1260000e-001
	  4.8000000e-007  1.3900000e-001
	  4.8500000e-007  1.6930000e-001
	  4.9000000e-007  2.0800000e-001
	  4.9500000e-007  2.5860000e-001
	  5.0000000e-007  3.2300000e-001
	  5.0500000e-007  4.0730000e-001
	  5.1000000e-007  5.0300000e-001
	  5.1500000e-007  6.0820000e-001
	  5.2000000e-007  7.1000000e-001
	  5.2500000e-007  7.9320000e-001
	  5.3000000e-007  8.6200000e-001
	  5.3500000e-007  9.1490000e-001
	  5.4000000e-007  9.5400000e-001
	  5.4500000e-007  9.8030000e-001
	  5.5000000e-007  9.9500000e-001
	  5.5500000e-007  1.0002000e+000
	  5.6000000e-007  9.9500000e-001
	  5.6500000e-007  9.7860000e-001
	  5.7000000e-007  9.5200000e-001
	  5.7500000e-007  9.1540000e-001
	  5.8000000e-007  8.7000000e-001
	  5.8500000e-007  8.1630000e-001
	  5.9000000e-007  7.5700000e-001
	  5.9500000e-007  6.9490000e-001
	  6.0000000e-007  6.3100000e-001
	  6.0500000e-007  5.6680000e-001
	  6.1000000e-007  5.0300000e-001
	  6.1500000e-007  4.4120000e-001
	  6.2000000e-007  3.8100000e-001
	  6.2500000e-007  3.2100000e-001
	  6.3000000e-007  2.6500000e-001
	  6.3500000e-007  2.1700000e-001
	  6.4000000e-007  1.7500000e-001
	  6.4500000e-007  1.3820000e-001
	  6.5000000e-007  1.0700000e-001
	  6.5500000e-007  8.1600000e-002
	  6.6000000e-007  6.1000000e-002
	  6.6500000e-007  4.4600000e-002
	  6.7000000e-007  3.2000000e-002
	  6.7500000e-007  2.3200000e-002
	  6.8000000e-007  1.7000000e-002
	  6.8500000e-007  1.1900000e-002
	  6.9000000e-007  8.2000000e-003
	  6.9500000e-007  5.7000000e-003
	  7.0000000e-007  4.1000000e-003
	  7.0500000e-007  2.9000000e-003
	  7.1000000e-007  2.1000000e-003
	  7.1500000e-007  1.5000000e-003
	  7.2000000e-007  1.0000000e-003
	  7.2500000e-007  7.0000000e-004
	  7.3000000e-007  5.0000000e-004
	  7.3500000e-007  4.0000000e-004
	  7.4000000e-007  3.0000000e-004
	  7.4500000e-007  2.0000000e-004
	  7.5000000e-007  1.0000000e-004
	  7.5500000e-007  1.0000000e-004
	  7.6000000e-007  1.0000000e-004
	  7.6500000e-007  0.0000000e+000
	  7.7000000e-007  0.0000000e+000];
  
  lu = interp1(tab(:,1), tab(:,2), lam/1E9);