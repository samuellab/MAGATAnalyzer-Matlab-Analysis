function x = randomWalkInLogarithms(n, sigma, minval, maxval)
%function randomWalkInLogarithms(n, sigma, minval, maxval)
%
%n is number of samples, minval is minimum value in log units(default: log(1) = 0), maxval is
%maximum value (default: log(255) = 5.5413)
%sigma is step size (in log units) (default log(255)/255 = 0.0217)

existsAndDefault('sigma', 0.0217);
existsAndDefault('minval', 0);
existsAndDefault('maxval', log(255));
%x0 = log(exp(maxval) + exp(minval));
x0 = 0.5 * (maxval + minval);


x = zeros([1 n]);
dx = sigma*randn([1 n]);
x(1) = x0;
for j = 2:length(dx)
    x(j) = x(j-1) + dx(j);
    if (x(j) < minval)
        x(j) = 2*minval-x(j);
    end
    if (x(j) > maxval)
        x(j) = 2*maxval - x(j);
    end
end

x = randomRound(exp(x));

function z = randomRound(y)
z = floor(y);
z = z+(rand(size(z)) < z-y);


% 
% 
% x = x0 + cumsum(sigma*randn([1 n]));
% ind = find (x > maxval | x < minval, 1, 'first');
% while (~isempty(ind))
%     if (x(ind) < minval)
%         x(ind:end) = 2*minval - x(ind:end);
%     else
%         x(ind:end) = 2*maxval - x(ind:end);
%     end
%     ind = find (x > maxval | x < minval, 1, 'first');
% end
% 
% x = round(exp(x));
