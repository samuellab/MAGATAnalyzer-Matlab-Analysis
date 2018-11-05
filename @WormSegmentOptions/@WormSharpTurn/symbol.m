function [sym,color] = symbol(st)
%function [sym,color] = symbol(st)
%gets a 1-letter latex code for the type symbol

if (length(st) == 1)
    if (~isfinite(st.userCode))
        switch(st.typeCode)
            case -1 
                sym = '$\Omega$';
                color = 'm';
            case 0
                sym = 'b';
                color = 'g';
            case {1,2}
               sym = 'R';
               color = 'r';
        end
    else
        switch(st.userCode)
            case -1 
                sym = '$\mathbf{\Omega$}';
                color = 'm';
            case 0
                sym = '$\mathbf{RR}$';
                color = 'g';
            case 1
               sym = '$\mathbf{R}$';
               color = 'r';
            case 2
                sym = '$\mathbf{b}$';
                color = 'g';
            case 100
                sym = '?';
                color = 'y';        
            case 250
                sym = '\o';
                color = 'c';
            case 254
                sym = 'mt';
                color = 'y';
        end
    end
    return;
end

for j = 1:length(st)
    [s,c] = sym(j).symbol;
    sym{j} = s;
    color{j} = c;
end
