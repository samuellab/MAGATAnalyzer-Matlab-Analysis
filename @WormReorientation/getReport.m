function  str = getReport(reo, varargin)
%generates a report about a run
%function  getReport(reo, varargin)
%
%output: STR (a string, or cell of strings if multiple runs)
%input: reo < WormReorientation

precision = '%.1f';
varargin = assignApplicable(varargin);

infostr = ['Reorienation: duration ' num2str(diff(reo.track.getDerivedQuantity('eti', false, [reo.startInd reo.endInd])), precision) ...
           's.  Consists of ' num2str(reo.numTurns) ' turns.'];
ststr = 'Turn sequence: ';
for j = 1:length(reo.sharpTurn)
    ststr = [ststr reo.sharpTurn(j).type ' '];
end

degstr = ['Theta in: ', num2str(rad2deg(reo.thetaIn), precision) ', Theta out: ' num2str(rad2deg(reo.thetaOut), precision) ', D Theta: ' ...
    num2str(rad2deg(reo.dTheta), precision)];

str = {infostr, ststr, degstr};
