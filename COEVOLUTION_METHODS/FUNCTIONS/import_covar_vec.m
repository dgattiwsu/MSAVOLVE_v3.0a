function import_covar_vec(filename)
% Import the file

rawData1 = importdata(filename);

% For some simple files (such as a CSV or JPEG files), IMPORTDATA might
% return a simple array.  If so, generate a structure so that the output
% matches that from the Import Wizard.

[~,name] = fileparts(filename);
newData1.(genvarname(name)) = rawData1;

% Create new variables in the base workspace from those fields. Remember to
% change the workspace to 'caller' if this function is used inside another
% function.

vars = fieldnames(newData1);
for i = 1:length(vars)
    assignin('base', vars{i}, newData1.(vars{i}));
end

