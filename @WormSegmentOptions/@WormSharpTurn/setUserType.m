function success = setUserType (st, charCode) 

if ~ischar(charCode)
    warningMessage();
    success = false;
    return;
end
success = true;
switch (upper(charCode(1)))
    case 'O' %mega turn
        st.userCode = -1;
    case 'D' %ouble reverse
        st.userCode = 0;
    case 'R' %eversal
        st.userCode = 1;
    case {'B','P'} %blip or pause
        st.userCode = 2;
    case 'C' %an't tell
        st.userCode = 100;
    case 'N' %ot a turn
        st.userCode = 250;
    case 'M' %ultiple turns 
        st.userCode = 254;
    otherwise
        warningMessage();
        success = false;
end

function warningMessage()

disp ('You did not enter a valid label for this turn.  Choices are:')
disp ('[O]mega turn - forms a ball');
disp ('[R]eversal - changes direction about 180 degrees without curling up');
disp ('[D]ouble reverse - changes direction twice within this single turn without curling up');
disp ('[M]ultiple turns - multiple turns (other than the double reverse) incorrectly grouped together here');
disp ('[P]ause - it just stops moving for a moment or two');
disp ('[C]an''t tell what''s going on here');
disp ('[N]ot a turn - this should not have been flagged as a turn at all');



        