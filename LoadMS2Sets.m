function Data=LoadMS2Sets(DataType)

%Loads all data sets of a certain type and outputs them into the structure
%Data

%DataType is the tab in the XLS file. The code figures out which XLS file
%and folders to use.


%Get some of the default folder
[SourcePath,FISHPath,DefaultDropboxFolder,MS2CodePath,SchnitzcellsFolder]=...
    DetermineLocalFolders;

%MS2Pausing folder:
if strcmp(DataType,'hbBAC')|strcmp(DataType,'Eve2')|strcmp(DataType,'snaBAC')|...
        strcmp(DataType,'snaBACNoPrimary')|strcmp(DataType,'snaBACNoShadow')|strcmp(DataType,'P2PPausing')|...
        strcmp(DataType,'hbNoPrimary')|strcmp(DataType,'kniBAC')|strcmp(DataType,'hbNoShadow')|...
        strcmp(DataType,'kniNoPrimary')|strcmp(DataType,'kniBAC')|strcmp(DataType,'kniNoShadow')|...
        strcmp(DataType,'snaNoPrimary')|strcmp(DataType,'snaNoShadow');
    [SourcePath,FISHPath,DropboxFolder,MS2CodePath,SchnitzcellsFolder]=...
        DetermineLocalFolders('2014-03-15-HbBACA');
    PausingXLSName='DataStatusPausing.xlsx';
elseif strcmp(DataType,'zld')|strcmp(DataType,'Shelby')
    [SourcePath,FISHPath,DropboxFolder,MS2CodePath,SchnitzcellsFolder]=...
        DetermineLocalFolders('2013-09-09-zld-X1');
    PausingXLSName='DataStatus.xlsx';
elseif strcmp(DataType,'MCP-GFP 5'' 2-Spot')
    [SourcePath,FISHPath,DropboxFolder,MS2CodePath,SchnitzcellsFolder]=...
        DetermineLocalFolders('2014-06-29-MCP(2)-X1(2S)')
    PausingXLSName='Data Status.xlsx';
elseif strcmp(DataType,'MCP-GFP 5'' Data')
    error('Take care of this')
else
    error('Add this data type to the code')
end



[StatusNum,StatusTxt]=xlsread([DropboxFolder,filesep,PausingXLSName],DataType);

CompileRow=find(strcmp(StatusTxt(:,1),'AnalyzeLiveData Compile Particles'));
CompiledSets=find(strcmp(StatusTxt(CompileRow,:),'READY')|strcmp(StatusTxt(CompileRow,:),'ApproveAll'));

clear SetNames
clear APDivisions
clear MeanFits
clear MeanFitsUp
clear Schnitzcells

for i=1:length(CompiledSets)
    SetName=StatusTxt{6,CompiledSets(i)};
    Quotes=strfind(SetName,'''');
    Prefix=SetName((Quotes(1)+1):(Quotes(end)-1));
    
    %Need to try this in case there's some incompatibility in terms of the
    %structures. This is because we might have data sets that have been
    %compiled using different versions of CompileParticles.m
    try
        Data(i)=load([DropboxFolder,filesep,Prefix,filesep,'CompiledParticles.mat']);
    catch
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
                warning(['Getting rid of field ', FieldsToCopy{j}])
            end
            Data(i)=DataTemp;
        else
            error('Need to take care of these other cases')
        end
    end
    
    if exist([DropboxFolder,filesep,Prefix,filesep,'APDivision.mat'])
        APDivisions(i)=load([DropboxFolder,filesep,Prefix,filesep,'APDivision.mat']);
    else
        warning('APDivisions.mat not found. This is a stupid way to check. Have the code check if this experiment is DV or AP instead')
    end
        
       
    %Fit results assuming the same slopes
    if exist([DropboxFolder,filesep,Prefix,filesep,'MeanFits.mat'])
        MeanFits(i)=load([DropboxFolder,filesep,Prefix,filesep,'MeanFits.mat']);
    else
        warning('MeanFits.mat not found. This is a stupid way to check. Have the code check if this experiment is DV or AP instead')
    end
    
    Schnitzcells(i)=load([DropboxFolder,filesep,Prefix(1:end),filesep,Prefix(1:end),'_lin.mat']);
    SetNames{i}=SetName;
    
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
    Ellipses(i)=load([DropboxFolder,filesep,Prefix,filesep,'Ellipses.mat']);
    %Count ellipses
    if exist([DropboxFolder,filesep,Prefix,filesep,'CountEllipses.mat'])
        CountEllipses(i)=load([DropboxFolder,filesep,Prefix,filesep,'CountEllipses.mat']);
    else
        %This is not the smarter way to do this. It relies on having at the
        %end a set that has been analyzed
        %CountEllipses(i)=[];
    end
end

%Pad the APDivisions in case we're missing it for some sets
if length(APDivisions)<length(CompiledSets)
    APDivisions(length(CompiledSets)).APDivision=[];
end


%Now add the SetName and APDivision information
for i=1:length(Data)
    Data(i).SetName=SetNames{i};
    
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
   
    Data(i).schnitzcells=Schnitzcells(i).schnitzcells;
    Data(i).Ellipses=Ellipses(i).Ellipses;
end

