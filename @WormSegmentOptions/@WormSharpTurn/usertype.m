function typeCode = usertype(st)
%function typeCode = usertype(st)
%
%gets the type strings based on the type codes for sharp turns

if (length(st) == 1)
    if (~isfinite(st.userCode))
        typeCode = 'not labeled';
        return;
    end
    switch(st.userCode)
        case -1 
            typeCode = 'omega turn';
        case 0
            typeCode = 'double reverse';
        case 1
            typeCode = 'reversal';
        case 2
            typeCode = 'blip or pause';
        case 100
            typeCode = 'could not tell';
        case 250
            typeCode = 'not a sharp turn'; 
        case 254
            typeCode = 'multiple turns';
    end
    return;
end

for j = 1:length(st)
    typeCode{j} = st(j).usertype;
end
