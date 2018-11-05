function [ Y ] = heavisideStep( X )
%HEAVISIDEFUNCTION Summary of this function goes here
%   We are using this so we can compile the code since MATLAB does not
%   allow compilation of functions of the Symbolic Math Toolbox

Y = zeros(size(X));
Y(X > 0) = 1;
Y(X == 0) = .5;

end

