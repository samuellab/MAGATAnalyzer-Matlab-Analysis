

[adim,po] = blank8x10Figure(fignum);
pp = [adim.lx3, adim.h0 - adim.h, adim.w3, adim.h];
poss = {pp, pp, pp};

poss{2}(1) = adim.cx3;
poss{3}(1) = 2*poss{2}(1) - poss{1}(1);
poss{3}(2) = poss{1}(2);

poss2 = poss; % 2nd row of figure (3 columns each)
poss2{1}(1) = poss{1}(1);
poss2{1}(2) = adim.cx3 + 0.15;
poss2{2}(2) = poss2{1}(2);
poss2{3}(2) = poss2{1}(2);

positions = {poss2{1}, poss2{2}};

j=1;


h{j}.pos = positions{j};
h{j}.ax = axes('Position', positions{j});
b = pcolor(taus, dts, -LogLBerlin.LL.logP_fitR + min(min(LogLBerlin.LL.logP_fitR)) ); shading interp; axis square
clear caxis
caxis(gca, [-100 0]);
colorbar
b.Parent.XAxis.TickValues = [0.5 10 20 30];
b.Parent.XAxis.TickLabels = {'0', '10', '20', '30'};

set(gcf, 'renderer', 'opengl')

j = 2;
h{j}.pos = positions{j};
h{j}.ax = axes('Position', positions{j});
a = pcolor(taus, dts, -LogL42a.LL.logP_fitR + min(min(LogL42a.LL.logP_fitR)) ); shading interp; axis square
a.Parent.XAxis.TickValues = [0.5 10 20 30];
a.Parent.XAxis.TickLabels = {'0', '10', '20', '30'};
clear caxis
caxis(gca, [-100 0]);
a = colorbar; a.Visible = 'off';

set(gcf, 'renderer', 'opengl')