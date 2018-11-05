function colorWheelHSV()
%function colorWheelHSV()
%
%make a color wheel

dh = 0.01;
h = 0:dh:1;
x = cos(h*2*pi);
x1 = cos((h+dh)*2*pi);
y = sin(h*2*pi);
y1 = sin((h+dh)*2*pi);

for j = 1:length(h)      
    patch([0 x(j) x1(j) 0], [0 y(j) y1(j) 0], hsv2rgb(h(j), 1,1), 'EdgeColor', 'none'); hold on;
end

axis equal;
