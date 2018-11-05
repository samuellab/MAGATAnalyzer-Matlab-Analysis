

exponential = 0;
new_axes = 0; % if you want to define new [dt,tau] grid and interpolate

conf_region = 1; % draw confidence ellipse, with radii found with fit2DLL.m


t=-pi:0.01:pi;

x0_42a=11; % tau_best
y0_42a=.7; % dt_best
tau=5; % tau confidence radius
dt=.15; % dt confidence radius
x1=x0_42a+tau*cos(t);
y1=y0_42a+dt*sin(t);

x0_ber=6; 
y0_ber=.85;
tau=2.5;
dt=.07;
t=-pi:0.01:pi;
xb1=x0_ber+tau*cos(t);
yb1=y0_ber+dt*sin(t);

LLs = [LogL_Berlin, LogL_42a];

clear dts taus
for k=1:2

LLstruct = LLs(k);
% LLstruct_fine = LLs_fine(k);

dts{k} = LLstruct.dts;
taus{k} = LLstruct.taus;

nll{k} = LLstruct.LL.logP_fitR;
zz{k} =  -nll{k} + min(min(nll{k}));


a1 = [.05 .1 .2 .3 .4];
a3 = [1.1:.4:1.9, 3];
a2 = .5:.2:1;
dts_new = [a1 a2 a3];
taus_new = taus{1};

[lxx, lxy] = meshgrid(taus{k}, dts{k});
[lxx2, lxy2] = meshgrid(taus_new, dts_new);
lxdata = [lxx(:) lxy(:)];

lxdata2 = [lxx2(:) lxy2(:)];
nll_new{k} = interp2(lxx, lxy, nll{k}, lxx2, lxy2, 'nearest');
zz_new{k} =  -nll_new{k} + min(min(nll_new{k}));
end

if(new_axes)
    zz_42a = zz_new{2};
    zz_berlin = zz_new{1};
    tauaxis_42a = taus_new;
    tauaxis_ber = taus_new;
    dtaxis_42a = dts_new;
    dtaxis_ber = dts_new;
else
    zz_42a = zz{2};
    zz_berlin = zz{1};
    tauaxis_42a = taus{2};
    dtaxis_42a = dts{2};
    tauaxis_ber = taus{1};
    dtaxis_ber = dts{1};
end
    
if(exponential)
    Z_42a = exp(zz_42a);
    Z_berlin = exp(zz_berlin);
    paxis = [0 1];
    pticks = [0 0.2 0.4 0.6 0.8 1];
else
    Z_42a = zz_42a;
    Z_berlin = zz_berlin;
    paxis = [-14 0];
    pticks = [-14 -7 0];
end


subplot(1,2,2)
b1 = pcolor(tauaxis_42a, dtaxis_42a, Z_42a ); shading interp; axis square
if(conf_region)
    hold on; plot(x1, y1, 'w:', 'LineWidth', 2);% plot(x2, y2, 'w:', 'LineWidth', 1.5); plot(x3, y3, 'w:', 'LineWidth', 1.5); plot(x4, y4, 'w:', 'LineWidth', 1.5);
    plot(x0_42a, y0_42a, 'w.', 'MarkerSize', 20);
    hold off;
end
b1.Parent.YAxis.Limits = [min(dtaxis_42a) 2];
b1.Parent.XAxis.Limits = [min(tauaxis_42a) 25];
b1.Parent.XAxis.TickValues = [0.5 10 20 30];
b1.Parent.XAxis.TickLabels = {'0', '10', '20', '30'};
clear caxis
caxis(gca, paxis);
a1 = colorbar; a1.Visible = 'off';
subplot(1,2,1)
b2 = pcolor(tauaxis_ber, dtaxis_ber, Z_berlin ); shading interp; axis square
if(conf_region)
    hold on; plot(xb1, yb1, 'w:', 'LineWidth', 2);
    plot(x0_ber, y0_ber, 'w.', 'MarkerSize', 20);
    hold off;
end
b2.Parent.YAxis.Limits = [min(dtaxis_ber) 2];
b2.Parent.XAxis.Limits = [min(tauaxis_ber) 20];
clear caxis
caxis(gca, paxis );
a2 = colorbar;
b2.Parent.XAxis.TickValues = [0.5 10 20 30];
b2.Parent.XAxis.TickLabels = {'0', '10', '20', '30'};
set(gcf, 'renderer', 'painters');

a2.Position(1) = .47;
a1.Position(1) = .9;
b1.Parent.Position(1) = .65;
b2.Parent.Position(1) = .06;
b1.Parent.Position(2:3) = b2.Parent.Position(2:3);
a2.Position(2) = .29;
a2.Position(4) = .45;
a2.Position(2) = .3;
a2.Position(4) = .43;
a2.Ticks = pticks;
a2.FontSize = 12;
b1.Parent.FontSize = 12;
b2.Parent.FontSize = 12;
b1.Parent.XAxis.TickValues = [0.5 10 20 30];
b2.Parent.XAxis.TickValues = [0.5 10 20 30];
b1.Parent.YAxis.TickValues = [0.5 1 1.5 2 2.5 3];
b2.Parent.YAxis.TickValues = [0.5 1 1.5 2 2.5 3];