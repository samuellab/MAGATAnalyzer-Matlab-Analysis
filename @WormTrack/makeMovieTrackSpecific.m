function mo = makeMovieTrackSpecific(track, mo, pt, iind) 

if (mo.labelSharpTurns && ~isempty(track.sharpTurn))
    xl = get(mo.Axes(1), 'XLim');
    yl = get(mo.Axes(1), 'YLim');
    r = [xl yl];
    st = track.sharpTurn;
    start = [st.startInd];
    stop = [st.endInd];
    st = st(stop >= iind - mo.ptbuffer & start <= iind + mo.ptbuffer);
    loc = double(st.getDerivedQuantity(mo.locField, 'position', 'centralInd'));
    validinds = (insideRect(r, loc));
    st = st(validinds);
    loc = loc(:,validinds);
    s = find([st.endInd] < iind, 1, 'last');
    if (isempty(s))
        s = 1;
    end
    e = find([st.startInd] > iind, 1, 'first');
    if (isempty(e))
        e = length(st);
    end
    s = max(s - mo.stbuffer, 1);
    e = min(e + mo.stbuffer, length(st));
    for j = s:e
        [s,c] = st(j).symbol;
        text(loc(1,j), loc(2,j),s, 'Color', c, 'Parent', mo.Axes(1), 'HorizontalAlignment', 'center', ...
            'VerticalAlignment', 'middle');
    end
    
    ind = find(st.containsIndex(iind), 1, 'first');
    if ~isempty(ind)
        [~,c] = st(ind).symbol;
        if (isfinite(st(ind).userCode))
            str = [st(ind).usertype ' (hand)'];
        else
            str = [st(ind).type ' (auto)'];
        end
        text(xl(1), yl(2), str, 'Color', c, 'Parent', mo.Axes(1), 'HorizontalAlignment', 'left', ...
            'VerticalAlignment', 'top');
    end
end
       
        