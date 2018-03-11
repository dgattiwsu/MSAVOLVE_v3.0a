function import_omes_matrix_caller(fileToRead1)
%  Imports data from the specified file

newData1 = importdata(fileToRead1);

% Create new variables in the caller workspace from those fields.

vars = fieldnames(newData1);
for i = 1:length(vars)
    assignin('caller', vars{i}, newData1.(vars{i}));
end

