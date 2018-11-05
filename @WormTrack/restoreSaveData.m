function track = restoreSaveData(track)
%function track = restoreSaveData(track)

sd = track.saveData;

 
if isempty(sd) || ~all(isfield(sd, {'sharpTurnStart', 'sharpTurnEnd', 'sharpTurnTypeCode', 'sharpTurnUserCode'})) || isempty(sd.sharpTurnStart) || isempty(sd.sharpTurnEnd) || isempty(sd.sharpTurnTypeCode)
    return;
end


for j = 1:length(sd.sharpTurnStart)
   sharpTurn(j) = WormSharpTurn(track, sd.sharpTurnStart(j), sd.sharpTurnEnd(j)); %#ok<*AGROW>
   sharpTurn(j).userCode = sd.sharpTurnUserCode(j);
   if (sign(sharpTurn(j).typeCode) ~= sign(sd.sharpTurnTypeCode(j)))
       disp (['on load, calculated sharp turn type code ' num2str(sharpTurn(j).typeCode) ...
           ' disagrees with saved type code ' num2str(sd.sharpTurnTypeCode(j))]);
   end
   sharpTurn(j).typeCode = sd.sharpTurnTypeCode(j);
   
end

track.sharpTurn = sharpTurn;
track.segmentTrack(track.so, 'UseExistingSharpTurns', true);
%{
sd.sharpTurnStart = [track.sharpTurn.startInd];
sd.sharpTurnEnd = [track.sharpTurn.endInd];
sd.sharpTurnTypeCode = [track.sharpTurn.typeCode];
sd.sharpTurnUserCode = [track.sharpTurn.userCode];
track.saveData = sd;
%}
