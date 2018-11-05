function ind = binarySearch(x, xval, le, ind1, ind2)
%DOES NOT WORK AS YOU MIGHT EXPECT - AMBIGUITIES IN HOW EDGES ARE HANDLEDE
%function ind = binarySearch(x, xval, le, ind1, ind2)
%finds xval in ascending sorted list x between ind1 and ind2
%designed to be used to find nonoverlapping range
%le = left edge
%if le is true x(ind) <= xval < x(ind+1)
%if le is false xval < x(ind+1);
if (ind2 - ind1 < 2)
    if (le)
        if (xval <= x(ind1))
            ind = ind1;
            return;
        else
            if (xval <= x(ind2))
                ind = ind2;
            else
                ind = [];
            end
            return;
        end
    else
        if (x(ind2) < xval)
            ind = ind2;
            return;
        else
            if (x(ind1) < xval)
                ind = ind1;
            else
                ind = [];
            end
            return;
        end
%         if (xval < x(ind1))
%             ind = ind1;
%             return;
%         else
%             if (xval < x(ind2))
%                 ind = ind2;
%             else
%                 ind = [];
%             end
%             return;
%         end
    end
end
    
ii = ceil((ind1 + ind2)/2);
if (le)
    if (x(ii) <= xval)
        ind = binarySearch(x, xval, le, ii, ind2);
        return;
    else
        ind = binarySearch(x, xval, le, ind1, ii);
        return;
    end
else
    if (x(ii) < xval)
        ind = binarySearch(x, xval, le, ii, ind2);
        return;
    else
        ind = binarySearch(x, xval, le, ind1, ii);
        return;
    end
end
