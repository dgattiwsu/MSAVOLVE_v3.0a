function import_coev_mat(filename)

% Import the file
rawData1 = importdata(filename);

% For some simple files (such as a CSV or JPEG files), IMPORTDATA might
% return a simple array.  If so, generate a structure so that the output
% matches that from the Import Wizard.

[~,name] = fileparts(filename);
newData1.(genvarname(name)) = rawData1;

% Create new variables in the base workspace from those fields. If called
% from inside a function, remember to replace the 'base' workspace with the
% 'caller' workspace.

vars = fieldnames(newData1);
for i = 1:length(vars)
    assignin('base', vars{i}, newData1.(vars{i}));
end

