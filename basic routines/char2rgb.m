function rgbvec = char2rgb (charcolor)
%function rgbvec = char2rgb (charcolor)
%
%converts a character color (one of 'r','g','b','c','m','y','k','w') to a 3
%value RGB vector
%if charcolor is a string (vector of chars), the result is a Nx3 matrix of
%color values, where N is the length of charcolor

if (~exist('charcolor','var') || ~ischar(charcolor))
    warning('RGB2VEC:NOTC', 'You must pass a character (rgbcmykw)');
    rgbvec = [];
    return;
end
if (~ischar(charcolor))
    rgbvec = charcolor;
    return;
end
rgbvec = zeros(length(charcolor), 3);
charwarning = false;
for j = 1:length(charcolor)
    switch(lower(charcolor(j)))
        case 'r'
            rgbvec(j,:) = [1 0 0];
        case 'g'
            rgbvec(j,:) = [0 1 0];
        case 'b'
            rgbvec(j,:) = [0 0 1];
        case 'c'
            rgbvec(j,:) = [0 1 1];
        case 'm'
            rgbvec(j,:) = [1 0 1];
        case 'y'
            rgbvec(j,:) = [1 1 0];
        case 'w'
            rgbvec(j,:) = [1 1 1];
        case 'k'
            rgbvec(j,:) = [0 0 0];
        otherwise
            charwarning = true;
    end
end

if (charwarning)
    warning('RGB2VEC:BADC', 'Only r,g,b,c,m,y,k,and w are recognized colors');
end
        