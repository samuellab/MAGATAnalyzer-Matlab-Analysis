function typeCode = type(st)
%function typeCode = type(st)
%
%gets a string interpretation of the typeCode

if (length(st) == 1)
    switch(st.typeCode)
        case -1 
            typeCode = 'omega turn';
        case 0
            typeCode = 'double reverse or blip';
        case 1
            typeCode = 'reversal';
        case 2
            typeCode = 'second reversal';
    end
    return;
end

for j = 1:length(st)
    typeCode{j} = st(j).type;
end
