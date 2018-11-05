function testAssignApplicapleAndFixTypes

boolval = false;
numval = 10;
charval = 'foobar';


assignApplicableAndFixTypes({'boolval', 'true'});
boolval
whos boolval

assignApplicableAndFixTypes({'boolval', '0'});
boolval
whos boolval

assignApplicableAndFixTypes({'numval', '147'});
numval
whos numval

assignApplicableAndFixTypes({'charval', 'true'});
charval
whos charval