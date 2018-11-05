function nas = fieldNamesAndSizes (structin)

fn = fieldnames(structin);

for j = 1:length(fn)
    bob = structin.(fn{j}); 
    rob = whos('bob');
    nas(j).name = fn{j};
    nas(j).bytes = rob.bytes;
end
[~,I] = sort([nas.bytes], 'descend');
nas = nas(I);
