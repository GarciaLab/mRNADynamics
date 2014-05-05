function Data=LoadMS2Sets(DataType)

%Loads all data sets of a certain type and outputs them into the structure
%Data

%DataType is the tab in the XLS file. The code figures out which XLS file
%and folders to use.


%Get some of the default folder
[SourcePath,FISHPath,DefaultDropboxFolder,MS2CodePath,SchnitzcellsFolder]=...
    DetermineLocalFolders;

%MS2Pausing folder:
if strcmp(DataType,'hbBAC')|strcmp(DataType,'Eve2')
    [SourcePath,FISHPath,DropboxFolderPausing,MS2CodePath,SchnitzcellsFolder]=...
        DetermineLocalFolders('2014-03-15-HbBACA');
    PausingXLSName='DataStatusPausing.xlsx';
else
    error('Add this data type to the code')
end



[StatusNum,StatusTxt]=xlsread([DropboxFolderPausing,filesep,PausingXLSName],DataType);

CompileRow=find(strcmp(StatusTxt(:,1),'AnalyzeLiveData Compile Particles'));
CompiledSets=find(strcmp(StatusTxt(CompileRow,:),'READY')|strcmp(StatusTxt(CompileRow,:),'ApproveAll'));

clear SetNames
clear APDivisions
clear MeanFits
clear Schnitzcells

for i=1:length(CompiledSets)
    SetName=StatusTxt{6,CompiledSets(i)};
    Quotes=strfind(SetName,'''');
    Prefix=SetName((Quotes(1)+1):(Quotes(end)-1));
    
    %Need to try this in case there's some incompatibility in terms of the
    %structures. This is because we might have data sets that have been
    %compiled using different versions of CompileParticles.m
    try
        Data(i)=load([DropboxFolderPausing,filesep,Prefix,filesep,'CompiledParticles.mat']);
    catch
        %If this fails figure out what's the missing field
        DataTemp=load([DropboxFolderPausing,filesep,Prefix,filesep,'CompiledParticles.mat']);
        
        FieldNamesData=fieldnames(Data);
        FieldNamesDataTemp=fieldnames(DataTemp);
        
        %If there are new fields we'll just get rid of them here
        if length(FieldNamesData)<length(FieldNamesDataTemp)
            %Figure out which fields to copy
            FieldsToCopy=FieldNamesDataTemp(~ismember(FieldNamesDataTemp,FieldNamesData));
            
            for j=1:length(FieldsToCopy)
                DataTemp=rmfield(DataTemp,FieldsToCopy{j});
                warning(['Getting rid of field ', FieldsToCopy{j}])
            end
            Data(i)=DataTemp;
        else
            error('Need to take care of these other cases')
        end
    end
        
    APDivisions(i)=load([DropboxFolderPausing,filesep,Prefix,filesep,'APDivision.mat']);
    %Fit results assuming the same slopes
    %MeanFits(i)=load([DropboxFolderPausing,filesep,Prefix,filesep,'MeanFits.mat']);
    Schnitzcells(i)=load([DropboxFolderPausing,filesep,Prefix(1:end),filesep,Prefix(1:end),'_lin.mat']);
    SetNames{i}=SetName;
    
    if exist([DropboxFolderPausing,filesep,Prefix,filesep,'FitIntegralResults.mat'])
        IntegralFits(i)=load([DropboxFolderPausing,filesep,Prefix,filesep,'FitIntegralResults.mat']);
    end
    
    if exist([DropboxFolderPausing,filesep,Prefix,filesep,'IndividualFits.mat'])
        IndividualFits(i)=load([DropboxFolderPausing,filesep,Prefix,filesep,'IndividualFits.mat']);
    end
    
    %Load Ellipses
    Ellipses(i)=load([DropboxFolderPausing,filesep,Prefix,filesep,'Ellipses.mat']);
    %Count ellipses
    if exist([DropboxFolderPausing,filesep,Prefix,filesep,'CountEllipses.mat'])
        CountEllipses(i)=load([DropboxFolderPausing,filesep,Prefix,filesep,'CountEllipses.mat']);
    else
        %This is not the smarter way to do this. It relies on having at the
        %end a set that has been analyzed
        %CountEllipses(i)=[];
    end
end

%Now add the SetName and APDivision information
for i=1:length(Data)
    Data(i).SetName=SetNames{i};
    Data(i).APDivision=APDivisions(i).APDivision;
    %Data(i).MeanFits=MeanFits(i).FitResults;
    
    if exist('IntegralFits')
        Data(i).IntegralFits=IntegralFits(i).FitResults;
    end
    
    if exist('IndividualFits')
        Data(i).IndividualFits=IndividualFits(i).FitResultsIndiv;
    end
    
   
    Data(i).schnitzcells=Schnitzcells(i).schnitzcells;
    Data(i).Ellipses=Ellipses(i).Ellipses;
end

