
function [mean, var] = BayesianVarEstimate(s0, Ddt, StimVal)
% function [mean, var] = BayesianVarEstimate(s0, Ddt, StimVal)
% Evaluate Bayesian-optimal estimate of variance, see DeWeese&Zador: https://www.mitpressjournals.org/doi/pdf/10.1162/089976698300017403

clear Ps meanvar;

dS = 0.5;
sigmas = 0.05:dS:5; %sigma axis for integration

[S1, S2] = meshgrid(sigmas, sigmas);

Nxs = @(x, s) normpdf(x, 0, s); %(3.22) from paper linked above - Distribution of stimuli given stddev; P(x|s)
Nss = @(s1,s2) normpdf(s1, s2, sqrt(2*Ddt)); % (3.21) - Diffusion prior on stddev P(sigma1 | sigma2), D=diffusion constant, dt = time step


%First time step after change from slow to shigh

Psx_0 = @(s) Nxs(StimVal(1), s).* Nss(s, s0);
omega = dS*trapz( Psx_0(sigmas) );
Ps(1,:) = Psx_0(sigmas)/omega;

%Psx = Psx_0;

Psx = Psx_0(sigmas)/omega;
% tic
for n=2:length(StimVal)
        
        Pss_n = dS*trapz(Nss(S1, S2).*repmat(Psx, length(sigmas), 1), 2);
        Psx_n = Nxs(StimVal(n), sigmas).*Pss_n';
        omega = dS*trapz(Psx_n);
        Psx = Psx_n/omega;
        Ps(n,:) = Psx;

end
% toc

for n=1:length(StimVal)
    mean(n) = trapz( sigmas.* Ps(n,:) ) / trapz(Ps(n,:));
    var(n) = trapz( sigmas.^2.* Ps(n,:) ) / trapz(Ps(n,:)) - mean(n)^2 ;
end

% shadedErrorPlot(1:n, meanvar, svar);

return




dX = 0.5;
stims = -10:dX:10; %x axis for integration
[X, S] = meshgrid(stims, sigmas);

if(xfirst)
    
    %First time step after change from slow to shigh
    
    Psx_0 = @(x, s) Nxs(x,s).*Nss(s, slow);
    omega = dS*trapz( Psx_0(X,S), 1);   
    Ps(1,:) = dX*trapz( Nxs(X, shigh).*Psx_0(X,S)./omega, 2);    
    Psx = Psx_0;
    
    tic
    for n=2:N
        
        % (2.3) - P(s_i | x_(j<i) )
        Pss_n = integral(@(a) Nss(S, a).*Psx(X, a)./omega, 0.05, 4.95, 'ArrayValued', true); 
        
        % (3.23) - P(s_i | x_(j<=i) )
        Psx_n = @(x, s) Nxs(x, s).* Pss_n ; 
        
        %normalization - Int( (3.23) * ds)
        omega = dS*trapz(Psx_n(X,S), 1); 
        
        % P(s_i) = Int( P(s|x) * P(x|s_high) * dx ) = average of P(s|x) over high-s stimulus distribution
        Ps(n,:) = dX*trapz(Nxs(X, shigh).*Psx_n(X,S)./omega, 2); 
        
        % update P(s | x)
        Psx = Psx_n; 
    end
    toc
    
    for n=1:size(Ps, 1)
        % <s> = Int( s * P(s) * ds )  = mean sigma
        meanvar(n) = trapz( sigmas.*Ps(n,:) )/trapz(Ps(n,:));
    end
    
    
    
    
    
    
else
   
    Psx_0 = @(x, s) Nxs(x,s).*Nss(s, slow);
    omega = dS*trapz(Psx_0(X,S), 1);
    Px(1,:) = dS*trapz(S.*Psx_0(X,S), 1);
    meanvar(1) = trapz(Nxs(stims, shigh).*Px(1,:)./omega)/trapz(Nxs(stims, shigh));
    Psx = Psx_0;
    
    tic
    for n=2:N
        
        % (2.3) - P(s_i | x_(j<i) )
        Pss_n = integral(@(a) Nss(S, a).*Psx(X, a)./omega, 0.05, 4.95, 'ArrayValued', true); 
        
        % (3.23) - P(s_i | x_(j<=i) )
        Psx_n = @(x, s) Nxs(x, s).* Pss_n ; 
        
        %normalization - Int( (3.23) * ds)
        omega = dS*trapz(Psx_n(X,S), 1); 
        
        % <s|x> = Int( P(s|x) * s * ds ) = mean sigma (s), given stim. (x)
        Px(n,:) = dS*trapz(S.*Psx_n(X,S), 1);
        
        %<s> = Int( <s|x> * P(x|s_high) * dx) = mean sigma
        meanvar(n) = trapz(Nxs(stims, shigh).*Px(n,:)./omega)/trapz(Nxs(stims, shigh));
        
        % update P(s | x)
        Psx = Psx_n; 
    end
    toc

    
end
    
plot(meanvar);
end
    