function [alpha, dt, sigAligned] = simpleLinearTimeRegistration(sig, sigRef)
%     sig = sig(:);
%     sigRef = sigRef(:);
    ncp = floor(.9 * length(sig));
    
    
    op = optimoptions('fmincon', 'Display', 'off');
    x = fmincon(@(x) -warpProduct(x(1), x(2)), [1 0], [], [], [], [], [.5 -length(sig)/4], [2 length(sig)/4],[],op);
   % x = fmincon(@(x) -warpProduct(x(1),0), [1], [], [], [], [], .5,2);    
    alpha = x(1);
    dt = x(2);
    %dt = 0;
    txf = 1:length(sig);
    tc1f = mean(txf) - dt;
    txw1f = 1/alpha * (txf - tc1f) + mean(tx);
    sigAligned = interp1(txw1f, sig, txf, 'linear', 0);
    
    function val = warpProduct(a, t)
        tx = 1:length(sig);
        tc1 = mean(tx) - t/2;
        tc2 = tc1 + t;
        txw1 = 1/sqrt(a) * (tx - tc1) + mean(tx);
        txw2 = sqrt(a) * (tx - tc2) + mean(tx);
        ti = mean(tx) - (ncp+1)/2 + (1:ncp);
        val = sum(sum(interp1(txw1, sig, ti, 'linear', 0).*interp1(txw2, sigRef, ti, 'linear', 0)));
    end
end