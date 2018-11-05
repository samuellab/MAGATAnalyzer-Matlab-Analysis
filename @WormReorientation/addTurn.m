function addTurn(reo, st)
%WormReorientation/addTurn
%adds a turn to the reorientation
reo.track = st.track;
reo.sharpTurn = [reo.sharpTurn st];
[~,I] = sort([reo.sharpTurn.startInd]);
reo.sharpTurn = reo.sharpTurn(I);

%reo.calculateMetrics();
