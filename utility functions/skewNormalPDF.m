function pdf = skewNormalPDF(x, u, s, a)
%function pdf = skewNormalPDF(x, u, s, a)
%skew-normal distribution with location parameter u, scale parameter s, and skew-parameter a
%http://en.wikipedia.org/wiki/Skew_normal_distribution
%pdf = 1/sqrt(2*pi)*exp(-(x-u).^2/(2*s^2)).*(1+erf(a * (x - u)/(sqrt(2)*s)));

%1/(pi*s)*exp(-(x-u)^2/(2*s^2)) * Sqrt[\[Pi]/2] Erfc[(a (u - x))/(Sqrt[2] s)]
%pdf = 1/(sqrt(2*pi).*s).*exp(-(x-u).^2./(2*s.^2)).*(1+erf(a .* (x - u)./((sqrt(2).*s))));


%pdf = 1./ sqrt(2* pi) ./ s .* exp(-(x-u).^2./(2*s.^2)) .* erfc(a .* (u - x)./sqrt(2)./s);

%Erfc[-((a*(-u + x))/(Sqrt[2]*s))]/(E^((-u + x)^2/(2*s^2))*Sqrt[2*Pi]*s)
pdf = erfc(-((a.*(-u + x))./(sqrt(2).*s))).*exp(-(-u + x).^2/(2.*s^2))./(sqrt(2*pi).*s);

end

