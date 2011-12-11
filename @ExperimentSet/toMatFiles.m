function toMatFiles(eset, fstub)
% saves every experiment in eset to a separate mat flile
% function toMatFiles(eset, fstub)
%
% saves every experiment in eset to a separate mat flile
% variable name is experiment_j in file fstub_experiment_j.mat
%
% outputs: none (to disk)
% inputs: 
%   ESET: a member of the ExperimentSet class
%   FSTUB: the file stub to which _experiment_1.mat, etc. will be appended

for j = 1:length(eset.expt)
    eval(['experiment_' num2str(j) ' = eset.expt(j);']);
    save ([fstub '_experiment_' num2str(j)], ['experiment_' num2str(j)]);
end
