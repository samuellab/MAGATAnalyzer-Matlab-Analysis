function a = surveyorFormula (pts)
%function a = surveyorFormula (pts)

a = 0;
for k = 1:(length(pts)-1)
    a = a + det([pts(:,k) pts(:,k+1)]);
end
a = a + det([pts(:,end) pts(:,1)]);
