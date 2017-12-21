function Data=LoadMS2Sets(DataType)
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
%No options
%
%OUTPUT
%Returns the Data structure containing all of the relevant datasets from your
%DataType tab in dataStatus.xlsx
%
%Author (contact): Hernan Garcia (hgarcia@berkeley.edu)
%Created: Unknown
%Last Updated: Unknown

%Get some of the default folders
[SourcePath,FISHPath,DefaultDropboxFolder,MS2CodePath,PreProcPath]=...
    DetermineLocalFolders;
[SourcePath,FISHPath,DropboxFolder,MS2CodePath,PreProcPath]=...
    DetermineLocalFolders;

%Now, get a list of all possible other Dropbox folders
[Dummy,XLS]=xlsread([MS2CodePath,filesep,'..',filesep,'ComputerFolders.XLSX']);
DropboxRows=find(~cellfun(@isempty,strfind(XLS(:,1),'Dropbox')));
DropboxFolders=XLS(DropboxRows,2);

%Look in DataStatus.XLSX in each DropboxFolder and find the tab given by
%the input variable DataType.
DataStatusToCheck=[];
for i=1:length(DropboxFolders)
    DDataStatus=dir([DropboxFolders{i},filesep,'DataStatus.*']);
    if length(DDataStatus)>1
        error(['More than one DataStatus.XLS found in folder ',DropboxFolders{i}])
    elseif length(DDataStatus)==1
        [Dummy,Sheets] = xlsfinfo([DropboxFolders{i},filesep,DDataStatus(1).name]);
        FindSheets=strcmpi(Sheets,lower(DataType));
        if sum(FindSheets)==1
            DataStatusToCheck=[DataStatusToCheck,i];
        end
    end
end
%Check that there aren't two DataStatus files with the same tab name
if length(DataStatusToCheck)>1
    error(['More than one DataStatus.XLSX found with a tab named ',DataType])
elseif length(DataStatusToCheck)==0
    error(['No DataStatus.XLSX found with a tab named ',DataType])
end

%Redefine the DropboxFolder according to the DataStatus.XLSX we'll use
DropboxFolder=DropboxFolders{DataStatusToCheck};
%Now, load the DataStatus.XLSX
D=dir([DropboxFolder,filesep,'DataStatus.*']);
[StatusNum,StatusTxt]=xlsread([DropboxFolder,filesep,D(1).name],DataType);


%Which data sets are approved?
CompileRow=find(strcmp(StatusTxt(:,1),'AnalyzeLiveData Compile Particles')|...
    strcmp(StatusTxt(:,1),'CompileParticles')|...
    strcmp(StatusTxt(:,1),'CompileNuclearProtein'));
CompiledSets=find(strcmp(StatusTxt(CompileRow,:),'READY')|strcmp(StatusTxt(CompileRow,:),'ApproveAll'));

if isempty(CompiledSets)
    error('No ApproveAll or READY sets found')
end

clear SetNames
clear APDivisions
clear MeanFits
clear MeanFitsUp
clear Schnitzcells



%Check the consistency between all the data acquired and analyzed in terms
%of ExperimentType, ExperimentAxis, and APResolution. Get this out of
%MovieDatabase.xlsx
[XLSNum,XLSTxt,XLSRaw]=xlsread([DefaultDropboxFolder,filesep,'MovieDatabase.xlsx']);
ExperimentTypeColumn=find(strcmp(XLSRaw(1,:),'ExperimentType'));
ExperimentAxisColumn=find(strcmp(XLSRaw(1,:),'ExperimentAxis'));
APResolutionColumn = find(strcmp(XLSRaw(1,:),'APResolution'));
DataFolderColumn=find(strcmp(XLSRaw(1,:),'DataFolder'));

ExperimentType=[];
ExperimentAxis=[];
APResolution=[];

PrefixRow=find(strcmp(StatusTxt(:,1),'Prefix:'));
for i=1:length(CompiledSets)
    SetName=StatusTxt{PrefixRow,CompiledSets(i)};
    Quotes=strfind(SetName,'''');
    Prefix=SetName((Quotes(1)+1):(Quotes(end)-1));

    Dashes=findstr(Prefix,'-');
    PrefixRowMovieDatabase=find(strcmp(XLSRaw(:,DataFolderColumn),[Prefix(1:Dashes(3)-1),'\',Prefix(Dashes(3)+1:end)]));
        if isempty(PrefixRowMovieDatabase)
            PrefixRowMovieDatabase=find(strcmp(XLSRaw(:,DataFolderColumn),[Prefix(1:Dashes(3)-1),'/',Prefix(Dashes(3)+1:end)]));
            if isempty(PrefixRowMovieDatabase)
                error('Could not find data set in MovieDatabase.XLSX. Check if it is defined there.')
            end
        end
    if isempty(PrefixRowMovieDatabase)
        error('Entry not found in MovieDatabase.xlsx')
    end
    
    %Load and check the experiment details consistency
    if ~isempty(ExperimentType)
        if ~strcmp(ExperimentType,XLSRaw{PrefixRowMovieDatabase,ExperimentTypeColumn})
           error('Inconsistent experiment types found among the data sets.') 
        end
    else
        ExperimentType=XLSRaw{PrefixRowMovieDatabase,ExperimentTypeColumn};
    end
    if ~isempty(ExperimentAxis)
        if ~strcmp(ExperimentAxis,XLSRaw{PrefixRowMovieDatabase,ExperimentAxisColumn})
           error('Inconsistent experiment axis found among the data sets.') 
        end
    else
        ExperimentAxis=XLSRaw{PrefixRowMovieDatabase,ExperimentAxisColumn};
    end
    if ~isempty(APResolution)
        if APResolution~=XLSRaw{PrefixRowMovieDatabase,APResolutionColumn}
           error('Inconsistent axis resolution found among the data sets.') 
        end
    else
        APResolution = XLSRaw{PrefixRowMovieDatabase,APResolutionColumn};
    end
    
    
    
end




%Find and load the different prefixes
PrefixRow=find(strcmp(StatusTxt(:,1),'Prefix:'));
for i=1:length(CompiledSets)
    SetName=StatusTxt{PrefixRow,CompiledSets(i)};
    SetNames{i}=SetName;
    Quotes=strfind(SetName,'''');
    Prefix=SetName((Quotes(1)+1):(Quotes(end)-1));
    
    %Load CompiledParticles if it exists. This constitutes the main part of
    %the Data output. However, we will later add more information to this
    %structure as well as use it to check the consistency of the analysis
    %performed with the different data sets.
    
    if exist([DropboxFolder,filesep,Prefix,filesep,'CompiledParticles.mat'])
        %Need to try this in case there's some incompatibility in terms of the
        %structures. This is because we might have xls sets that have been
        %compiled using different versions of CompileParticles.m
        try
            DataTemp=load([DropboxFolder,filesep,Prefix,filesep,'CompiledParticles.mat']);
            DataTemp=orderfields(DataTemp);
            Data(i)=DataTemp;
        catch
            warning('Incompatible analysis in the data sets. Re-check this part of the code.')
            %If this fails figure out what's the missing field
            DataTemp=load([DropboxFolder,filesep,Prefix,filesep,'CompiledParticles.mat']);
            FieldNamesData=fieldnames(Data);
            FieldNamesDataTemp=fieldnames(DataTemp);

            %If there are new fields we'll just get rid of them here
            if length(FieldNamesData)<length(FieldNamesDataTemp)
                %Figure out which fields to copy
                FieldsToCopy=FieldNamesDataTemp(~ismember(FieldNamesDataTemp,FieldNamesData));

                for j=1:length(FieldsToCopy)
                    DataTemp=rmfield(DataTemp,FieldsToCopy{j});
                    DataTemp=orderfields(DataTemp);
                    warning(['Getting rid of field ', FieldsToCopy{j}])
                end
                Data(i)=DataTemp;
            else
                FieldsToCopy = FieldNamesData(~ismember(FieldNamesData,FieldNamesDataTemp));
                for j=1:length(FieldsToCopy)
                    DataTemp.(FieldsToCopy{j}) = {};               
                    DataTemp=orderfields(DataTemp);
                end
            end
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
        
        
        
        if exist([DropboxFolder,filesep,Prefix,filesep,'APDivision.mat'])
            APDivisions(i)=load([DropboxFolder,filesep,Prefix,filesep,'APDivision.mat']);
        else
            warning('APDivisions.mat not found. This is a stupid way to check. Have the code check if this experiment is DV or AP instead')
        end


        %Fit results assuming the same slopes
        if exist([DropboxFolder,filesep,Prefix,filesep,'MeanFits.mat'],'file')
            MeanFits(i)=load([DropboxFolder,filesep,Prefix,filesep,'MeanFits.mat']);
        else
            %Commented out by AR 12/11/17 since this comment seems
            %outdated.
%             warning(['MeanFits.mat not found for ',Prefix,...
%                 '. This is a stupid way to check. Have the code check if this experiment is DV or AP instead.'])
        end

        try
            Schnitzcells(i)=load([DropboxFolder,filesep,Prefix(1:end),filesep,Prefix(1:end),'_lin.mat']);
        catch
            warning('_lin.mat not found.');
        end
        

        %Fit to the integrals
        if exist([DropboxFolder,filesep,Prefix,filesep,'FitIntegralResults.mat'])
            IntegralFits(i)=load([DropboxFolder,filesep,Prefix,filesep,'FitIntegralResults.mat']);
        end


        %Fits to individual traces
        if exist([DropboxFolder,filesep,Prefix,filesep,'IndividualFits.mat'])
            IndividualFits(i)=load([DropboxFolder,filesep,Prefix,filesep,'IndividualFits.mat']);
        end

        %Integrated amount accounting from degradation. This is generated using
        %Jacques' code
    %     if exist([DropboxFolder,filesep,Prefix,filesep,'AccumulationData.mat'])
    %         AccumulationData(i)=load([DropboxFolder,filesep,Prefix,filesep,'AccumulationData.mat']);
    %     end

        if exist([DropboxFolder,filesep,Prefix,filesep,'MeanFitsUp.mat'])
            MeanFitsUp(i)=load([DropboxFolder,filesep,Prefix,filesep,'MeanFitsUp.mat']);
        end


        %Load Ellipses
        try
        Ellipses(i)=load([DropboxFolder,filesep,Prefix,filesep,'Ellipses.mat']);
        catch
            warning('Ellipses.mat not found.')
        end

        % Load Particles
        load([DropboxFolder,filesep,Prefix,filesep,'Particles.mat']);
        ParticleTemp(i).Particles=Particles;
        %Count ellipses
        if exist([DropboxFolder,filesep,Prefix,filesep,'CountEllipses.mat'])
            CountEllipses(i)=load([DropboxFolder,filesep,Prefix,filesep,'CountEllipses.mat']);
        else
            %This is not the smarter way to do this. It relies on having at the
            %end a set that has been analyzed
            %CountEllipses(i)=[];
        end
    end
    
    
    %Load CompiledNuclei if it exists
    if exist([DropboxFolder,filesep,Prefix,filesep,'CompiledNuclei.mat'])
        DataNuclei(i)=load([DropboxFolder,filesep,Prefix,filesep,'CompiledNuclei.mat']);
    end
end


%Now add the SetName and APDivision information
if exist([DropboxFolder,filesep,Prefix,filesep,'CompiledParticles.mat'])
    for i=1:length(Data)
        Data(i).SetName=SetNames{i};

        if exist('ImageRotation')
            Data(i).ImageRotation=ImageRotation(i);
        end
        
        if exist('APDivisions')
            Data(i).APDivision=APDivisions(i).APDivision;
        end

        if exist('IntegralFits')
            if i<=length(IntegralFits)
                Data(i).IntegralFits=IntegralFits(i).FitResults;
            end
        end

        if exist('IndividualFits')
            if i<=length(IndividualFits)
                Data(i).IndividualFits=IndividualFits(i).FitResultsIndiv;
            end
        end

        if exist('AccumulationData')
            if i<=length(AccumulationData)
                Data(i).AccumulationData=AccumulationData(i).AcumData;
            end
        end

        if exist('MeanFits')
            if i<=length(MeanFits)
                Data(i).MeanFits=MeanFits(i).FitResults;
            end
        end

       if exist('MeanFitsUp')
            if i<=length(MeanFitsUp)
                Data(i).MeanFitsUp=MeanFitsUp(i).FitResults;
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
if exist('DataNuclei')
    for i=1:length(DataNuclei)
        DataNuclei(i).SetName=SetNames{i};
    end
end



%If we have both particles and nuclei, then combine everything
if exist('Data')&exist('DataNuclei')
    DataTemp=Data;
    clear Data
    for i=1:length(DataTemp)
        Data(i).Particles=DataTemp(i);
        Data(i).Nuclei=DataNuclei(i);
    end
elseif (~exist('Data'))&exist('DataNuclei')
    Data=DataNuclei;
elseif  (~exist('Data'))&(~exist('DataNuclei'))
    error('No CompiledParticles found. Check DynamicsResults folder as well as DataStatus.XLSX.')
end
    


