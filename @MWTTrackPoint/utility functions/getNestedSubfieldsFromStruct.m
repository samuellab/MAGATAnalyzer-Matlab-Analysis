function val = getNestedSubfieldsFromStruct (st, dotnamedfield)
%function val = getNestedSubfieldsFromStruct (st, dotnamedfield)
%what st.(dotnamedfield) would give you if that were legal
%ex: dotnamefield = 'superstruct.substruct'
    subfields = regexp(dotnamedfield, '\.', 'split');
    ds = [st.(subfields{1})];
    for j = 2:length(subfields)
        ds = [ds.(subfields{j})];
    end
    val = ds;
end