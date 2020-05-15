function [Data, prefixes, resultsFolder,...
    ignoredPrefixes, dataTypeTabContents] = LoadMS2Sets(dataType, varargin)
%
% Data = LoadMS2Sets(DataType)
%
% DESCRIPTION
% Loads all data sets of a certain type and outputs them into the structure
% Data
%
% PARAMETERS
% DataType: This is a string that is identical to the name of the tab in
% dataStatus.xlsx that you wish to analyze.
%
% OPTIONS
% 'noCompiledNuclei'
% 'justPrefixes'
% 'localMovieDatabase' - use MovieDatabase in same local directory as
% DataStatus file
% 'dataStatusFolder', dataStatusFolder: 
%
% OUTPUT
% Data: Returns the Data structure containing all of the relevant datasets from your
%           DataType tab in dataStatus.xlsx
% prefixes: List of prefixes associated with the tab
% resultsFolder: location of the data
% ignoredPrefixes: list of prefixes not used to compile data
%


%Author (contact): Hernan Garcia (hgarcia@berkeley.edu)
%Created:
%Last Updated: 5/13/2020 MT (previously 3/18/2020 JL)

% Initialize
prefixes = {};
Data = struct();
resultsFolder = '';
dataStatusFilename = 'DataStatus.*';    %NB: This naming convention is now enforced inside findDataStatus.m

% Process the options
[noCompiledNuclei, justPrefixes, inputOutputFits, inputOutputModel, ... 
    localMovieDatabase,dataStatusFolder] = determineLoadMS2SetsOptions(varargin);

%Get some of the default folders
[rawDataPath, ProcPath, ~, MS2CodePath, PreProcPath,...
    configValues, movieDatabasePath, movieDatabaseFolder, movieDatabase, ...
    allDropboxFolders] =  DetermineLocalFolders;

% Find all DataStatus.xlsx files in all DropboxFolders
dataStatusFolders = findDataStatus(allDropboxFolders);

%Look in all DataStatus.xlsx files to find the tab specified by the dataType
%user input
dataStatusWithDataTypeFolder = findDataStatusTab(dataStatusFolders, dataType);

%Redefine the DropboxFolder according to the DataStatus.xlsx we'll use
dropboxFolder = dataStatusWithDataTypeFolder;
resultsFolder = dropboxFolder;

%Create a subfolder in the resutlsFolder for the project specified by the
%dataType variable
projectFolder = [resultsFolder, filesep, dataType];
if ~exist(projectFolder,'dir')
    try
        mkdir(projectFolder);
    catch
        warning(['Failed to create new directory ', projectFolder])
    end
end

%Redefine the MovieDatabase according to the DropboxFolder we're using if
%using the LocalMovieDatabase option
if localMovieDatabase
    movieDatabasePath = [dropboxFolder,'\MovieDatabase.csv'];
    movieDatabase = csv2cell(movieDatabasePath, 'fromfile');
end


%Now, load the contents of the DataStatus.XLSX tab we just found
dataStatusDir = dir([dropboxFolder,filesep,dataStatusFilename]);
dataTypeTabContents = readcell([dropboxFolder,filesep,dataStatusDir(1).name],'Sheet', dataType);


%Which data sets are approved?
CompileRow=find(strcmpi(dataTypeTabContents(:,1),'AnalyzeLiveData Compile Particles')|...
    strcmpi(dataTypeTabContents(:,1),'CompileParticles')|...
    strcmpi(dataTypeTabContents(:,1),'CompileNuclearProtein'));
CompiledSets=find(strcmpi(dataTypeTabContents(CompileRow,:),'READY')|strcmpi(dataTypeTabContents(CompileRow,:),'ApproveAll'));

ignoredSets=find( (~strcmpi(dataTypeTabContents(CompileRow,:),'READY')) & (~strcmpi(dataTypeTabContents(CompileRow,:),'ApproveAll')));
ignoredSets = ignoredSets(2:end);

PrefixRow=strcmp(dataTypeTabContents(:,1),'Prefix:');

ignoredPrefixes = cell(1, length(ignoredSets));

for i=1:length(ignoredSets)
    SetName=dataTypeTabContents{PrefixRow,ignoredSets(i)};
    Quotes=strfind(SetName,'''');
    ignoredPrefixes{i} = SetName((Quotes(1)+1):(Quotes(end)-1));
end
    
    
if isempty(CompiledSets)
    error('No ApproveAll or READY sets found')
end

clear SetNames
clear APDivisions
clear MeanFits
clear MeanLinearFits
clear MeanFitsUp
clear MeanLinearFitsUp
clear Schnitzcells
clear MeanFitsMCMC
clear InputOutputFits
clear SingleParticleFitsMCMC



%Check the consistency between all the data acquired and analyzed in terms
%of ExperimentType, ExperimentAxis, and APResolution. Get this out of
%MovieDatabase
ExperimentType=[];
ExperimentAxis=[];
APResolution=[];

prefixes = cell(1, length(CompiledSets));

for i=1:length(CompiledSets)
    SetName=dataTypeTabContents{PrefixRow,CompiledSets(i)};
    Quotes=strfind(SetName,'''');
    Prefix=SetName((Quotes(1)+1):(Quotes(end)-1));
    prefixes{i} = Prefix;
    
    [~, ExperimentTypeFromDatabase, ExperimentAxisFromDatabase, ~, ~, APResolutionFromDatabase, ~,...
        ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~] = getExperimentDataFromMovieDatabase(Prefix, movieDatabase);
    
    
    %Load and check the experiment details consistency
    if ~isempty(ExperimentType)
        if ~strcmpi(ExperimentType, ExperimentTypeFromDatabase)
            error('Inconsistent experiment types found among the data sets.')
        end
    else
        ExperimentType = ExperimentTypeFromDatabase;
    end
    if ~isempty(ExperimentAxis)
        if ~strcmpi(ExperimentAxis, ExperimentAxisFromDatabase)
            error('Inconsistent experiment axis found among the data sets.')
        end
    else
        ExperimentAxis = ExperimentAxisFromDatabase;
    end
    if ~isempty(APResolution)
        if APResolution ~= APResolutionFromDatabase
            error('Inconsistent axis resolution found among the data sets.')
        end
    else
        APResolution = APResolutionFromDatabase;
    end
end

%7/15/19 JL: DropboxFolder was already defined above, why do we need this
%line here? I'm running into issues where I'm trying to load data from a
%different folder than my normal DynamicsResults, and using
%readMovieDatabase defaults my DropboxFolder to my normal DynamicsResults,
%resulting in a failure to load the datasets I want.
%[~,~,DropboxFolder] = readMovieDatabase(Prefix, optionalResults);

%Find and load the different prefixes
PrefixRow=find(strcmpi(dataTypeTabContents(:,1),'Prefix:'));
for i=1:length(CompiledSets)
    
    
    SetName=dataTypeTabContents{PrefixRow,CompiledSets(i)};
    SetNames{i}=SetName;
    Quotes=strfind(SetName,'''');
    Prefix=SetName((Quotes(1)+1):(Quotes(end)-1));
    
    if ~justPrefixes
        
        %Load CompiledParticles if it exists. This constitutes the main part of
        %the Data output. However, we will later add more information to this
        %structure as well as use it to check the consistency of the analysis
        %performed with the different data sets.
        
        if exist([dropboxFolder,filesep,Prefix,filesep,'CompiledParticles.mat'],'file')
            %Need to try this in case there's some incompatibility in terms of the
            %structures. This is because we might have xls sets that have been
            %compiled using different versions of CompileParticles.m
            try
                DataTemp=load([dropboxFolder,filesep,Prefix,filesep,'CompiledParticles.mat']);
                DataTemp=orderfields(DataTemp);
                Data(i)=DataTemp;
            catch
                [Data, DataTemp] = addFields(Data, DataTemp);
                Data(i) = DataTemp;
            end
            
            
            %         %Load information about the rotation of the zoomed-in image.
            %         if exist([DropboxFolder,filesep,Prefix,filesep,'APDetection.mat'])
            %             APDetection=load([DropboxFolder,filesep,Prefix,filesep,'APDetection.mat']);
            %             try
            %                 ImageRotation(i)=APDetection.ImageRotation;
            %             catch
            %                 error('Image rotation information not found in APDetection.mat. Rerun AddParticlePosition.m')
            %             end
            %         elseif strcmpi(ExperimentAxis,'ap')|strcmpi(ExperimentAxis,'dv')
            %             error(['APDetection.mat not found despite this experiment being on the ',ExperimentAxis,' axis'])
            %         end

            
            %Fit results assuming the same slopes
            if exist([dropboxFolder,filesep,Prefix,filesep,'MeanFits.mat'],'file')
                MeanFits(i)=load([dropboxFolder,filesep,Prefix,filesep,'MeanFits.mat']);
            else
            end
            
            %Linear slope fit results
            if exist([dropboxFolder,filesep,Prefix,filesep,'MeanLinearFits.mat'],'file')
                MeanLinearFits(i)=load([dropboxFolder,filesep,Prefix,filesep,'MeanLinearFits.mat']);
            end
            
            % Fit results from the MeanFitAPAsymmetric.m
            if exist([dropboxFolder,filesep,Prefix,filesep,'MeanFitsV2.mat'],'file')
                MeanFitsV2(i)=load([dropboxFolder,filesep,Prefix,filesep,'MeanFitsV2.mat']);
            else
            end
            
            % Fit results from the MeanFitAPAsymmetric.m
            if exist([dropboxFolder,filesep,Prefix,filesep,'MeanFitsAsymmetric.mat'],'file')
                MeanFitsAsymmetric(i)=load([dropboxFolder,filesep,Prefix,filesep,'MeanFitsAsymmetric.mat']);
            else
            end
            
            
            % Fit results from the FitTiltedTrapezoids_4Clicks.m
            if exist([dropboxFolder,filesep,Prefix,filesep,'MeanFitsV3.mat'],'file')
                MeanFitsV3(i)=load([dropboxFolder,filesep,Prefix,filesep,'MeanFitsV3.mat']);
            else
            end
            
            %Fit results using MeanFitsMCMC
            if exist([dropboxFolder,filesep,Prefix,filesep,'MeanFitsMCMC.mat'],'file')
                MeanFitsMCMC(i)=load([dropboxFolder,filesep,Prefix,filesep,'MeanFitsMCMC.mat']);
            else
            end
            %Fit results using InputOutputFits
            if inputOutputFits
                if exist([dropboxFolder,filesep,Prefix,filesep,'InputOutputFits_',...
                        inputOutputModel,'.mat'],'file')
                    InputOutputFits(i)=load([dropboxFolder,filesep,Prefix,filesep,...
                        'InputOutputFits_',inputOutputModel,'.mat']);
                else
                end
            end
            
            %Single particle results using MeanFitsMCMC
            if exist([dropboxFolder,filesep,Prefix,filesep,'SingleParticleFitsMCMC.mat'],'file')
                SingleParticleFitsMCMC(i)=load([dropboxFolder,filesep,Prefix,filesep,'SingleParticleFitsMCMC.mat']);
            else
            end
            
            try
                Schnitzcells(i)=load([dropboxFolder,filesep,Prefix(1:end),filesep,Prefix(1:end),'_lin.mat'], 'schnitzcells');
            catch
            end
            
            
            %Fit to the integrals
            if exist([dropboxFolder,filesep,Prefix,filesep,'FitIntegralResults.mat'],'file')
                IntegralFits(i)=load([dropboxFolder,filesep,Prefix,filesep,'FitIntegralResults.mat']);
            end
            
            
            %Fits to individual traces
            if exist([dropboxFolder,filesep,Prefix,filesep,'IndividualFits.mat'],'file')
                IndividualFits(i)=load([dropboxFolder,filesep,Prefix,filesep,'IndividualFits.mat']);
            end
            
            %Integrated amount accounting from degradation. This is generated using
            %Jacques' code
            %     if exist([DropboxFolder,filesep,Prefix,filesep,'AccumulationData.mat'])
            %         AccumulationData(i)=load([DropboxFolder,filesep,Prefix,filesep,'AccumulationData.mat']);
            %     end
            
            if exist([dropboxFolder,filesep,Prefix,filesep,'MeanFitsUp.mat'],'file')
                MeanFitsUp(i)=load([dropboxFolder,filesep,Prefix,filesep,'MeanFitsUp.mat']);
            end
            
            if exist([dropboxFolder,filesep,Prefix,filesep,'MeanLinearFitsUp.mat'],'file')
                MeanLinearFitsUp(i)=load([dropboxFolder,filesep,Prefix,filesep,'MeanLinearFitsUp.mat']);
            end
            
            
            %Load Ellipses
            try
                Ellipses(i)=load([dropboxFolder,filesep,Prefix,filesep,'Ellipses.mat'], 'Ellipses');
            catch
            end
            
            % Load Particles
            load([dropboxFolder,filesep,Prefix,filesep,'Particles.mat']);
            ParticleTemp(i).Particles=Particles;
            %Count ellipses
            if exist([dropboxFolder,filesep,Prefix,filesep,'CountEllipses.mat'],'file')
                CountEllipses(i)=load([dropboxFolder,filesep,Prefix,filesep,'CountEllipses.mat']);
            else
                %This is not the smarter way to do this. It relies on having at the
                %end a set that has been analyzed
                %CountEllipses(i)=[];
            end
        end
        
        % Here are the fields that we need to load in case there's no
        % CompiledParticles. i.e. Only CompiledNuclei.mat exists.
        % Load APDivisions if it exists
        if exist([dropboxFolder,filesep,Prefix,filesep,'APDivision.mat'],'file')
            APDivisions(i)=load([dropboxFolder,filesep,Prefix,filesep,'APDivision.mat']);
        else
        end
        
        %Load CompiledNuclei if it exists
        if exist([dropboxFolder,filesep,Prefix,filesep,'CompiledNuclei.mat'],'file') & ~noCompiledNuclei
            try
                DataNucleiTemp=load([dropboxFolder,filesep,Prefix,filesep,'CompiledNuclei.mat']);
                DataNucleiTemp=orderfields(DataNucleiTemp);
                DataNuclei(i)=DataNucleiTemp;
            catch
                [DataNuclei, DataNucleiTemp] = addFields(DataNuclei, DataNucleiTemp);
                DataNuclei(i) = DataNucleiTemp;
            end
        end
        
    end
end

if ~justPrefixes
    %Now add the SetName and APDivision information
    if exist([dropboxFolder,filesep,Prefix,filesep,'CompiledParticles.mat'],'file')
        for i=1:length(Data)
            Data(i).SetName=SetNames{i};
            
            if exist('ImageRotation','var')
                Data(i).ImageRotation=ImageRotation(i);
            end
            
            try
                if exist('APDivisions','var')
                    Data(i).APDivision=APDivisions(i).APDivision;
                end
            end
            
            if exist('IntegralFits','var')
                if i<=length(IntegralFits)
                    Data(i).IntegralFits=IntegralFits(i).FitResults;
                end
            end
            
            if exist('IndividualFits','var')
                if i<=length(IndividualFits)
                    Data(i).IndividualFits=IndividualFits(i).FitResultsIndiv;
                end
            end
            
            if exist('AccumulationData','var')
                if i<=length(AccumulationData)
                    Data(i).AccumulationData=AccumulationData(i).AcumData;
                end
            end
            
            if exist('MeanFits','var')
                if i<=length(MeanFits)
                    Data(i).MeanFits=MeanFits(i).FitResults;
                end
            end
            
            if exist('MeanFitsV2','var')
                if i<=length(MeanFitsV2)
                    Data(i).MeanFitsV2=MeanFitsV2(i).FitResults;
                end
            end
            
            if exist('MeanFitsAsymmetric','var')
                if i<=length(MeanFitsAsymmetric)
                    Data(i).MeanFitsAsymmetric=MeanFitsAsymmetric(i).FitResults;
                end
            end
            
            if exist('MeanLinearFits','var')
                if i<=length(MeanLinearFits)
                    Data(i).MeanLinearFits=MeanLinearFits(i).FitResults;
                end
            end
            
            if exist('MeanLinearFitsUp','var')
                if i<=length(MeanLinearFitsUp)
                    Data(i).MeanLinearFitsUp=MeanLinearFitsUp(i).FitResults;
                end
            end
            if exist('MeanFitsMCMC','var')
                if i<=length(MeanFitsMCMC)
                    Data(i).MeanFitsMCMC=MeanFitsMCMC(i).MCMCresults;
                end
            end
            if exist('InputOutputFits','var')
                if i<=length(InputOutputFits)
                    Data(i).InputOutputFits=InputOutputFits(i).MCMCresults;
                end
            end
            if exist('SingleParticleFitsMCMC','var')
                if i<=length(SingleParticleFitsMCMC)
                    Data(i).SingleParticleFitsMCMC=SingleParticleFitsMCMC(i).MCMCresults;
                end
            end
            try
                Data(i).schnitzcells=Schnitzcells(i).schnitzcells;
                Data(i).Ellipses=Ellipses(i).Ellipses;
            catch
            end
            
            Data(i).Particles=ParticleTemp(i).Particles;
        end
    end
    
% Add information to DataNuclei if it exists
if exist('DataNuclei','var')
    for i=1:length(DataNuclei)
        DataNuclei(i).SetName=SetNames{i};
        try
            if exist('APDivisions','var')
                DataNuclei(i).APDivision=APDivisions(i).APDivision;
            end
        end
    end
end



%If we have both particles and nuclei, then combine everything
%First, check to see if we have empty structures and clear them if we do
if exist('Data','var') && isempty(fieldnames(Data))
    clear Data
end
if exist('DataNuclei','var') && isempty(fieldnames(DataNuclei))
    clear DataNuclei
end

if exist('Data','var') && exist('DataNuclei','var')
    DataTemp=Data;
    clear Data
    for i=1:length(DataTemp)
        Data(i).Particles=DataTemp(i);
        Data(i).Nuclei=DataNuclei(i);
    end
elseif (~exist('Data', 'var')) && exist('DataNuclei', 'var')
    Data=DataNuclei;
elseif  (~exist('Data','var')) && (~exist('DataNuclei','var'))
    error('No CompiledParticles found. Check DynamicsResults folder as well as DataStatus.XLSX.')
end

if noCompiledNuclei
    DataTemp = Data;
    clear Data
    for i=1:length(DataTemp)
        Data(i).Particles=DataTemp(i);
    end
end

end

save([projectFolder, filesep, 'prefixes.mat'], 'prefixes');


end