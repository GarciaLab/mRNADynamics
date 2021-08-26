% Convenience function to load all contents of a .mat file to workspace, while also
% a) listing names of loaded variables and b) returning a matfile handle for variable manipulation
function F = loadMatFile(filepath, loadWholeFile)
	if nargin < 2
		loadWholeFile = false;
	end
	fprintf('Loading %s file into workspace...\n', filepath);
	F = matfile(filepath);
	varlist = whos(F);
	if loadWholeFile
		command = sprintf('load(''%s'')', filepath);
		evalin('base', command);
		
		varnames = extractfield(varlist, 'name');
		varnamesStrArray = convertCharsToStrings(varnames);
		varnamesStr = sprintf('%s, ', varnamesStrArray);
		varnamesStr = varnamesStr(1:end-2); % remove trailing comma and space.

		fprintf('Loaded vars: %s.\n', varnamesStr);
	else
		fprintf('Not loading any vars, just returning matfile handle.\n');
	end

end
