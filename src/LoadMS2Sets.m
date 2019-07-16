function [Data, Prefixes, resultsFolder] = LoadMS2Sets(DataType, varargin)
%
%Data = LoadMS2Sets(DataType)
%
%DESCRIPTION
%Loads all data sets of a certain type and outputs them into the structure
%Data
%
%PARAMETERS
%DataType: This is a string that is identical to the name of the tab in
%dataStatus.xlsx that you wish to analyze.
%
%OPTIONS
%'dontCompare': Skip comparing experiment settings
%
%OUTPUT
%Returns the Data structure containing all of the relevant datasets from your
%DataType tab in dataStatus.xlsx
%
%Author (contact): Hernan Garcia (hgarcia@berkeley.edu)
%Created:
%Last Updated: 1/13/2018. AR

Prefixes = {};

optionalResults = '';
compareSettings = true;

for i= 1:length(varargin)
    if strcmpi(varargin{i},'optionalResults')
        optionalResults = varargin{i+1};
    end
    if strcmpi(varargin{i},'dontCompare')
        compareSettings = false;
    end
end

%Get some of the default folders
[rawDataPath, ProcPath, DefaultDropboxFolder, MS2CodePath, PreProcPath,...
    configValues, movieDatabasePath, movieDatabaseFolder, movieDatabase] =  DetermineLocalFolders;

%Now, get a list of all possible other Dropbox folders
DropboxRows=contains(configValues(:,1),'Dropbox');
DropboxFolders=configValues(DropboxRows,2);

%Look in DataStatus.XLSX in each DropboxFolder and find the tab given by
%the input variable DataType.
DataStatusToCheck=[];
for i=1:length(DropboxFolders)
    DDataStatus=dir([DropboxFolders{i},filesep,'DataStatus.*']);
    if length(DDataStatus)>1
        error(['More than one DataStatus.XLS found in folder ',DropboxFolders{i}])
    elseif length(DDataStatus)==1
        [~,Sheets] = xlsfinfo([DropboxFolders{i},filesep,DDataStatus(1).name]);
        FindSheets=strcmpi(Sheets,DataType);
        if sum(FindSheets)==1
            DataStatusToCheck=[DataStatusToCheck,i];
        end
    end
end
%Check that there aren't two DataStatus files with the same tab name
if length(DataStatusToCheck)>1
    error(['More than one DataStatus.XLSX found with a tab named ',DataType])
elseif isempty(DataStatusToCheck)
    error(['No DataStatus.XLSX found with a tab named ',DataType])
end

%Redefine the DropboxFolder according to the DataStatus.XLSX we'll use
DropboxFolder=DropboxFolders{DataStatusToCheck};
resultsFolder = DropboxFolder;

%Now, load the DataStatus.XLSX
D=dir([DropboxFolder,filesep,'DataStatus.*']);
[~,StatusTxt]=xlsread([DropboxFolder,filesep,D(1).name],DataType);


%Which data sets are approved?
CompileRow=find(strcmpi(StatusTxt(:,1),'AnalyzeLiveData Compile Particles')|...
    strcmpi(StatusTxt(:,1),'CompileParticles')|...
    strcmpi(StatusTxt(:,1),'CompileNuclearProtein'));
CompiledSets=find(strcmpi(StatusTxt(CompileRow,:),'READY')|strcmpi(StatusTxt(CompileRow,:),'ApproveAll'));

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



%Check the consistency between all the data acquired and analyzed in terms
%of ExperimentType, ExperimentAxis, and APResolution. Get this out of
%MovieDatabase
ExperimentType=[];
ExperimentAxis=[];
APResolution=[];

PrefixRow=strcmp(StatusTxt(:,1),'Prefix:');
for i=1:length(CompiledSets)
    SetName=StatusTxt{PrefixRow,CompiledSets(i)};
    Quotes=strfind(SetName,'''');
    Prefix=SetName((Quotes(1)+1):(Quotes(end)-1));
    
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
PrefixRow=find(strcmpi(StatusTxt(:,1),'Prefix:'));
for i=1:length(CompiledSets)
    
    
    SetName=StatusTxt{PrefixRow,CompiledSets(i)};
    SetNames{i}=SetName;
    Quotes=strfind(SetName,'''');
    Prefix=SetName((Quotes(1)+1):(Quotes(end)-1));
    Prefixes = [Prefixes, Prefix];
    
    %Load CompiledParticles if it exists. This constitutes the main part of
    %the Data output. However, we will later add more information to this
    %structure as well as use it to check the consistency of the analysis
    %performed with the different data sets.
    
    if exist([DropboxFolder,filesep,Prefix,filesep,'CompiledParticles.mat'],'file')
        %Need to try this in case there's some incompatibility in terms of the
        %structures. This is because we might have xls sets that have been
        %compiled using different versions of CompileParticles.m
        try
            DataTemp=load([DropboxFolder,filesep,Prefix,filesep,'CompiledParticles.mat']);
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
        
        
        
        if exist([DropboxFolder,filesep,Prefix,filesep,'APDivision.mat'],'file')
            APDivisions(i)=load([DropboxFolder,filesep,Prefix,filesep,'APDivision.mat']);
        else
            warning('APDivisions.mat not found.')
        end
        
        
        %Fit results assuming the same slopes
        if exist([DropboxFolder,filesep,Prefix,filesep,'MeanFits.mat'],'file')
            MeanFits(i)=load([DropboxFolder,filesep,Prefix,filesep,'MeanFits.mat']);
        else
            warning('MeanFits.mat not found');
        end
        
        %Linear slope fit results
        if exist([DropboxFolder,filesep,Prefix,filesep,'MeanLinearFits.mat'],'file')
            MeanLinearFits(i)=load([DropboxFolder,filesep,Prefix,filesep,'MeanLinearFits.mat']);
        end
        
        % Fit results from the MeanFitAPAsymmetric.m
        if exist([DropboxFolder,filesep,Prefix,filesep,'MeanFitsV2.mat'],'file')
            MeanFitsV2(i)=load([DropboxFolder,filesep,Prefix,filesep,'MeanFitsV2.mat']);
        else
            warning('MeanFitsV2.mat not found');
        end
        
        % Fit results from the FitTiltedTrapezoids_4Clicks.m
        if exist([DropboxFolder,filesep,Prefix,filesep,'MeanFitsV3.mat'],'file')
            MeanFitsV3(i)=load([DropboxFolder,filesep,Prefix,filesep,'MeanFitsV3.mat']);
        else
            warning('MeanFitsV3.mat not found');
        end
        
        %Fit results assuming the same slopes
        if exist([DropboxFolder,filesep,Prefix,filesep,'MeanFitsMCMC.mat'],'file')
            MeanFitsMCMC(i)=load([DropboxFolder,filesep,Prefix,filesep,'MeanFitsMCMC.mat']);
        else
            warning('MeanFitsMCMC.mat not found');
        end
        
        try
            Schnitzcells(i)=load([DropboxFolder,filesep,Prefix(1:end),filesep,Prefix(1:end),'_lin.mat'], 'schnitzcells');
        catch
            warning('_lin.mat not found.');
        end
        
        
        %Fit to the integrals
        if exist([DropboxFolder,filesep,Prefix,filesep,'FitIntegralResults.mat'],'file')
            IntegralFits(i)=load([DropboxFolder,filesep,Prefix,filesep,'FitIntegralResults.mat']);
        end
        
        
        %Fits to individual traces
        if exist([DropboxFolder,filesep,Prefix,filesep,'IndividualFits.mat'],'file')
            IndividualFits(i)=load([DropboxFolder,filesep,Prefix,filesep,'IndividualFits.mat']);
        end
        
        %Integrated amount accounting from degradation. This is generated using
        %Jacques' code
        %     if exist([DropboxFolder,filesep,Prefix,filesep,'AccumulationData.mat'])
        %         AccumulationData(i)=load([DropboxFolder,filesep,Prefix,filesep,'AccumulationData.mat']);
        %     end
        
        if exist([DropboxFolder,filesep,Prefix,filesep,'MeanFitsUp.mat'],'file')
            MeanFitsUp(i)=load([DropboxFolder,filesep,Prefix,filesep,'MeanFitsUp.mat']);
        end
        
        if exist([DropboxFolder,filesep,Prefix,filesep,'MeanLinearFitsUp.mat'],'file')
            MeanLinearFitsUp(i)=load([DropboxFolder,filesep,Prefix,filesep,'MeanLinearFitsUp.mat']);
        end
        
        
        %Load Ellipses
        try
            Ellipses(i)=load([DropboxFolder,filesep,Prefix,filesep,'Ellipses.mat'], 'Ellipses');
        catch
            warning('Ellipses.mat not found.')
        end
        
        % Load Particles
        load([DropboxFolder,filesep,Prefix,filesep,'Particles.mat']);
        ParticleTemp(i).Particles=Particles;
        %Count ellipses
        if exist([DropboxFolder,filesep,Prefix,filesep,'CountEllipses.mat'],'file')
            CountEllipses(i)=load([DropboxFolder,filesep,Prefix,filesep,'CountEllipses.mat']);
        else
            %This is not the smarter way to do this. It relies on having at the
            %end a set that has been analyzed
            %CountEllipses(i)=[];
        end
    end
    
    
    %Load CompiledNuclei if it exists
    if exist([DropboxFolder,filesep,Prefix,filesep,'CompiledNuclei.mat'],'file')
        DataNuclei(i)=load([DropboxFolder,filesep,Prefix,filesep,'CompiledNuclei.mat']);
        
        if exist([DropboxFolder,filesep,Prefix,filesep,'APDivision.mat'],'file')
            APDivisions(i)=load([DropboxFolder,filesep,Prefix,filesep,'APDivision.mat'], 'APDivision');
        else
            warning('APDivision.mat not found.')
        end
        
    end
end


%Now add the SetName and APDivision information
if exist([DropboxFolder,filesep,Prefix,filesep,'CompiledParticles.mat'],'file')
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
        try
            Data(i).schnitzcells=Schnitzcells(i).schnitzcells;
            Data(i).Ellipses=Ellipses(i).Ellipses;
        catch
            warning('No schnitzcells or ellipses');
        end
        
        Data(i).Particles=ParticleTemp(i).Particles;
    end
end

% Add information to DataNuclei if it exists
if exist('DataNuclei','var')
    for i=1:length(DataNuclei)
        DataNuclei(i).SetName=SetNames{i};
        
        if exist('APDivisions','var')
            DataNuclei(i).APDivision=APDivisions(i).APDivision;
        end
    end
end



%If we have both particles and nuclei, then combine everything
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

end