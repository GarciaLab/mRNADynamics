function [A, B, C] = readImarisEllipsoidTables(imarisStatisticsFolder)

	ellipsoidFilePrefix = 'Ellipsoid_Axis_Length_';

	function ellipsoidTable = readEllipsoidFile(suffix)
		ellipsesFilePath = [imarisStatisticsFolder, filesep, ellipsoidFilePrefix, suffix, '.csv'];
	    opts = detectImportOptions(ellipsesFilePath);
	    opts.VariableNamesLine = 3;

	    ellipsoidTable = readtable(ellipsesFilePath, opts);
        ellipsoidTable = rmmissing(ellipsoidTable);
        ellipsoidTable = removevars(ellipsoidTable, {'Unit', 'Category'});
	end

	A = readEllipsoidFile('A');
	B = readEllipsoidFile('B');
	C = readEllipsoidFile('C');
 end
