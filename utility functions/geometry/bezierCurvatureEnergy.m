function [cv2, grad_cv2, cv] = bezierCurvatureEnergy(cpts, s)
%function [cv2, grad_cv2] = bezierCurvatureEnergy(cpts, s)
%right now, grad_cv2 is only correct for quadratic bezier (3 cpts)
%because I am lazy right now, also only produces correct values for 
%uniformly spaced s from 0 to 1
%
[~,Mp] = bezierMatrix(size(cpts,2)-1, s);

du = cpts*Mp;
ddu = diff(du,[],2)*length(s); 
ddu = 0.5*(ddu(:,[1 1:end]) + ddu(:,[1:end end]));
v2 = sqrt(sum(du.^2));

cv = (du(1,:).*ddu(2,:) - du(2,:).*ddu(1,:))./(v2.^(1.5));

cv2 = sum(cv.^2)/length(s);

num = (mean(du(1,:).*ddu(2,:) - du(2,:).*ddu(1,:)))^2;
coeff = 6/length(s)*num./(v2.^4);
gradx = sum((ones(size(Mp,1),1)*(coeff.*(cpts(1,:)*Mp)).*Mp),2);
grady = sum((ones(size(Mp,1),1)*(coeff.*(cpts(2,:)*Mp)).*Mp),2);

grad_cv2 = [gradx';grady'];