function stat = test_ReadYaml()
% this function tests reading in the yaml file

stat.ok = 1;
stat.desc = '';
try
    stat.test_ReadYaml_SimpleStructure = test_ReadYaml_SimpleStructure();   
    stat.test_ReadYaml_DateTime = test_ReadYaml_DateTime();  
catch    
    stat.ok = 0;
    stat.desc  = 'Program crash';
end

end


function stat = test_ReadYaml_SimpleStructure()

stat.ok = 1;
stat.desc = '';
try
    s = ReadYaml('simple.yaml');
    
    ages = [s.age];
    
    if not(isequal([33 27], ages))  || not(all(ismember({'John Smith', 'Mary Smith'}, {s.name}) ))
        stat.desc  = ' Wrong values loaded';
        stat.ok = 0;
    end
    
catch   
       stat.desc  = 'Program crash';
       stat.ok = 0;
end


end

function stat = test_ReadYaml_DateTime()

stat.ok = 1;
stat.desc = '';
try
    s = ReadYaml('time.yaml');
    
    if ~isa(s.Data.B1_S_SW{1},'DateTime')
        stat.desc  = ' Wrong data type of datetick';
        stat.ok = 0;
    end
    if isa(s.Data.B1_S_SW{2},'DateTime')
        stat.desc  = ' Wrong data type of datetick';
        stat.ok = 0;
    end
catch
       stat.desc  = 'Program crash';
       stat.ok = 0;
end
end