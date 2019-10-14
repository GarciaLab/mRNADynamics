function compareNuclearFluo(tabNames,prefix)

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
[~,~,userDynResPath,~,~] = ...
    DetermineLocalFolders(prefix);

% 1. Cycle through the length of tabNames
ultimateCompiliationFluo = zeros(1,length(tabNames));
for i = 1:length(tabNames)
    % 2. With each tab, check which prefixes should be included. This would
    % be indicated with a 1 in the corresponding cell
    % Note for the future: Should the indicator be a yes instead of 1?
    project = tabNames{i};
    prefix_cell = compilePrefixNames(project,userDynResPath);
    numOfDataSets = length(prefix_cell);
    compiledAvgNuclearFluo = zeros(numOfDataSets);
    compiledTotalNumNuclei = zeros(numOfDataSets);
    for j = 1:numOfDataSets
        prefix = prefix_cell{j};
        [avgNuclearFluo, totalNumNuclei] = extractNuclearFluoDosage(prefix);
        compiledAvgNuclearFluo(j) = avgNuclearFluo;
        compiledTotalNumNuclei(j) = totalNumNuclei;
    end
    ultimateCompiliationFluo(i) = sum(compiledAvgNuclearFluo)/...
        sum(compiledTotalNumNuclei);
end

plot(1:length(tabNames),ultimateCompiliationFluo,'.','MarkerSize',20)
xticks(1:length(tabNames))
xticklabels(tabNames)
set(gca,'TickLabelInterpreter','none',...
    'XTickLabelRotation',45)%,...
    %'Yscale','log') 
end