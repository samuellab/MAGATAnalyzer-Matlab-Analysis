function track = fillOutSaveData(track)
%function fillOutSaveData(track)



if (~isempty(track.sharpTurn))
    try
        sd.sharpTurnStart = [track.sharpTurn.startInd];
        sd.sharpTurnEnd = [track.sharpTurn.endInd];
        sd.sharpTurnTypeCode = [track.sharpTurn.typeCode];
        sd.sharpTurnUserCode = [track.sharpTurn.userCode];
    catch me
        disp(me.getReport);
        sd.sharpTurnStart = [];
        sd.sharpTurnEnd = [];
        sd.sharpTurnTypeCode = [];
        sd.sharpTurnUserCode = [];
    end
else
    sd.sharpTurnStart = [];
    sd.sharpTurnEnd = [];
    sd.sharpTurnTypeCode = [];
    sd.sharpTurnUserCode = [];
end

track.saveData = sd;
