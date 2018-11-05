function c = pathToCellOfStrings(p)
%function c = pathToCellOfStrings(p)
%
%p = 'c:\foo\bar\fubar.foo'
%c =  'c:\'    'foo'    'bar'    'fubar.foo'
j = 1;
while (j < 100)
    [p, nm, ext] = fileparts(p);
    if (~isempty(nm))
        c{j} = [nm ext];
    else
        c{j} = p;
        break;
    end
    j = j+1;
end

c = c(end:-1:1);

