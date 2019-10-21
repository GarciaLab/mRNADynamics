function ultimateCompiliationFluo = ...
    compareNuclearFluo(tabNames,prefix,varargin)
close all;
% This function will be given a cell of tab names in DataStatus and will
% compare the values of the nuclear fluo of those constructs.
% Flow:
% 1. Cycle through the length of tabNames through each tab
% 2. With each tab, check if the row entry "Ran ExportDataForFISH" has a 1 in
% the current data set of the tab.
% 3. With the movies identified, use extractNuclearFluo to get the average
% and the total number of nuclei used.
% 4. Plot the result on a plot that is avg nuclear fluo vs construct.
% 5. Think about what the error bar should be...
k = 1;
displayFigures = 0;
while k <= length(varargin)
    if strcmp(varargin(k),'displayFigures')
        displayFigures = 1;
    elseif strcmp(varargin(k),'histogram')
        displayFigures = 'histogram';
    end
    k = k + 1;
end


[~,~,userDynResPath,~,~] = ...
    DetermineLocalFolders(prefix);

%note for the future:
% Some how find a way to figure out from preprocessed data or meta data if
% there are multiple embryos imaged. 

% 1. Cycle through the length of tabNames
numOfProjects = length(tabNames);
ultimateCompiliationFluo = zeros(1,length(tabNames));
for i = 1:numOfProjects
    % 2. With each tab, check which prefixes should be included. This would
    % be indicated with a 1 in the corresponding cell
    % Note for the future: Should the indicator be a yes instead of 1?
    project = tabNames{i};
    prefix_cell = compilePrefixNames(project,userDynResPath);
    disp(prefix_cell)
    numOfDataSets = length(prefix_cell);
    compiledSumNuclearFluo = zeros(numOfDataSets);
    compiledTotalNumNuclei = zeros(numOfDataSets);
    for j = 1:numOfDataSets
        prefix = prefix_cell{j};
        [sumNuclearFluo, totalNumNuclei] = ...
            extractNuclearFluoDosage(prefix,displayFigures);
        compiledSumNuclearFluo(j) = sumNuclearFluo;
        compiledTotalNumNuclei(j) = totalNumNuclei;
    end
    ultimateCompiliationFluo(i) = sum(compiledSumNuclearFluo)/...
        sum(compiledTotalNumNuclei);
end

plot(1:numOfProjects,ultimateCompiliationFluo,'.','MarkerSize',20)
xticks(1:numOfProjects)
xticklabels(tabNames)
xlim([0 numOfProjects+1])
set(gca,'TickLabelInterpreter','none',...
    'XTickLabelRotation',45);%,...
    %'Yscale','log') 
end