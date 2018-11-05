function st = firstTurn(reo) 
% get the first turn of a reorientation or vector of reorientations
% function st = firstTurn(reo) 
%
% outputs: 
%   ST < sharpTurn;
% inputs:
%   REO < WormReorientation;

for j = 1:length(reo)
    st(j) = reo(j).sharpTurn(1); %#ok<AGROW>
end