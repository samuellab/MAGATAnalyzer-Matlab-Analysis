function testBezier

[X,Y] = getpts();
cpts = [X Y]';
pts = bezierFromControlPoints(cpts, 100);
plot (X,Y,'r.', X,Y,'k--',pts(1,:), pts(2,:), 'b-');
