% function rename_allvars_bystring(search_string,new_string). This script
% changes the name of all the variables of the workspace that contain a
% certain 'search_string' with new variables that contain a
% 'replace_string'. Unfortunately it cannot be used as a function because
% 'who' and following operations refer only to the function workspace and
% not to the 'global' workspace.

search_string = 'old_string';
replace_string = 'new_string';

allVars = who('-regexp',search_string,'global');

for i = 1:length(allVars)
newvarname = regexprep(allVars{i},search_string,replace_string);
eval( [ newvarname , ' = ', allVars{i},';']);
eval(['clear ', allVars{i},';']);
end

