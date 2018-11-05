function dkl = divergenceKL_perez_cruz_1D (px, qx)
%function dkl = divergenceKL_perez_cruz_1D (px, qx)
%
%estimates KL divergence between p & q, based on samples px and qx from p&q
%respectively
%
%algorithm from:1. F. Pérez-Cruz, in Information Theory, 2008. ISIT 2008. IEEE International Symposium on (IEEE, 2008; http://ieeexplore.ieee.org/xpls/abs_all.jsp?arnumber=4595271), pp. 1666–1670.
%does not work for me right now
%possible problem: finding epsilon < min (diff)
maxN = 1E6;
%mindx = (max(px)-min(px))/min(length(px)*10, 1E6); %keep axis size reasonable

pxx = unique(px)';
pxx = [pxx pxx(end)*1.01-pxx(end-1)*.01];
hh = histc(px, pxx);
pc = cumsum(hh);

paxis = linspace(pxx(1), pxx(end), maxN);
pe = interp1(pxx, pc/pc(end), paxis, 'linear');
dpe = [diff(pe) 0];

qxx = unique(qx)';

qxx = [qxx qxx(end)*1.01-qxx(end-1)*.01];
%qaxis = linspace(pxx(1), pxx(end), maxN);
hh = histc(qx, qxx);
qc = cumsum(hh);
qe = interp1(qxx, qc/qc(end), paxis, 'linear');
dqe = [diff(qe) 0];

plot (paxis, pe, paxis, qe)

dkl = mean(log(interp1(paxis, dpe, px)./interp1(paxis, dqe, px)));
