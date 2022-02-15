function [prefixes,prefixCellText] = getPrefixesFromDataStatusTab(dataTypeTabContents)

% function prefixes = getPrefixesFromDataStatusTab(dataTypeTabContents)
%
% DESCRIPTION
% Returns the Prefixes for all experiments in a project tab of 
% DataStatus.xlsx
%
% INPUT
% dataTypeTabContents: Cell array containing the contents of a
%                      DataStatus.xlsx tab
%
% OPTIONS
% N/A
%
% OUTPUT
% prefixes: n x 1 cell array containing the Prefixes for this project tab
% prefixCellText: n x 1 cell array of the full contents of the prefix cells
%                 from the DataStatus tab, only returning because of
%                 dependencies in LoadMS2Sets.m
%
% Author (contact): Meghan Turner (meghan_turner@berkeley.edu)
% Created: 5/17/2020
% Origin: Functionalized from code originally included in LoadMS2Sets.m,
%         written by Hernan Garcia (hggarcia@berkeley.edu)
% Last Updated: N/A

prefixes = {}; %#ok<*NASGU>
prefixCellText = {};

% Locate the prefix cells within the tab
prefixRow = find(strcmpi(dataTypeTabContents(:,1),'prefix:'));
prefixColumns = find( contains(...
    strrep( dataTypeTabContents(prefixRow,:) , ' ', ''), 'prefix', 'IgnoreCase', true) ); 

%the first column is a header
prefixColumns = prefixColumns(2:end);

if isempty(prefixRow) || isempty(prefixColumns)
    error('Unable to find Prefixes in DataStatus.xslx tab. Could be missing or incorrectly formatted.')
end
                 
prefixCellText = dataTypeTabContents(prefixRow,prefixColumns);

%Extract just the prefixes from the larger cell contents. The whole cell 
%will contain something like Prefix = 'YYYY-MM-DD-experimentName' but we 
%only want what's inside the single quotation marks
prefixes = cell(length(prefixColumns),1);
for i = 1:length(prefixCellText)
    currPrefixText = prefixCellText{i};
    prefixInQuotes = strfind(currPrefixText,'''');
    prefixName = currPrefixText((prefixInQuotes(1)+1):(prefixInQuotes(end)-1));
    prefixes{i} = prefixName;
end
