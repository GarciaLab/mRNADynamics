function StartTr2d(Prefix)

%Grab the folder associated with this Prefix
[RawDynamicsData,ProcessedData,DynamicsResults,MS2CodePath,PreProcessedData]=...
    DetermineLocalFolders(Prefix)

%Create a folder for the tr2d project inside PreProcessedData
mkdir([PreProcessedData,filesep,Prefix,filesep,'tr2dProject'])
%Create a folder inside the tr2d one for the exported segmentation and
%tracking
mkdir([PreProcessedData,filesep,Prefix,filesep,'tr2dProject',filesep,...
    'mRNADynamicsExport'])


%We need to convert the TIF sequence into a TIFF stack
%Make sure a TIF stack doesn't exist already
CreateStack=1;
if exist([PreProcessedData,filesep,Prefix,filesep,'tr2dProject',...
        filesep,'RAW.tiff'])
    Answer=input('Nuclear data has already been exported for tr2d. Do you want to export again? (Y/n)','s');
    if ~strcmpi(Answer,'n')
        CreateStack=0;
    end
end

if CreateStack
    %Load all nuclear images
    D=dir([PreProcessedData,filesep,Prefix,filesep,'*His*.tif']);
    for i=1:length(D)
        Image=imread([PreProcessedData,filesep,Prefix,filesep,...
            D(i).name]);
        imwrite(Image,[PreProcessedData,filesep,Prefix,filesep,'tr2dProject',...
            filesep,'RAW.tif'],'WriteMode','append')
    end
end

%Pull out the pixel size (in um) and the list of time stamps (in seconds).
%Save it as a CSV file that will serve a input for tr2d.
load([DynamicsResults,filesep,Prefix,filesep,'FrameInfo.mat'])
PixelSize=FrameInfo(1).PixelSize;
TimeStamps=[FrameInfo.Time];
%Save the information as CSV
CSVOutput(2,:)=TimeStamps;
CSVOutput(1,1)=PixelSize;
csvwrite([PreProcessedData,filesep,Prefix,filesep,'tr2dProject',...
    filesep,'mRNADynamics.csv'],CSVOutput);


%%%%%%%
%RUN TR2D
%%%%%%%%

%Check that we have the tr2d results and import them
if exist([PreProcessedData,filesep,Prefix,filesep,'tr2dProject',filesep,...
    'mRNADynamicsExport',filesep,'tr2d_objects.csv'])&...
    exist([PreProcessedData,filesep,Prefix,filesep,'tr2dProject',filesep,...
    'mRNADynamicsExport',filesep,'tr2d_tracks.csv'])

    ImportTr2d(Prefix)
else
    error('No tr2d results found, did you export?')    
end








