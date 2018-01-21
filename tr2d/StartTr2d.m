function StartTr2d(Prefix)

%Grab the folder associated with this Prefix
[RawDynamicsData,ProcessedData,DynamicsResults,MS2CodePath,PreProcessedData]=...
    DetermineLocalFolders(Prefix);

%Create a folder for the tr2d project inside PreProcessedData
mkdir([PreProcessedData,filesep,Prefix,filesep,'tr2dProject']);
%Create a folder inside the tr2d one for the exported segmentation and
%tracking
mkdir([PreProcessedData,filesep,Prefix,filesep,'tr2dProject',filesep,...
    'mRNADynamicsExport']);


%We need to convert the TIF sequence into a TIFF stack
%Make sure a TIF stack doesn't exist already
CreateStack=1;
filenameRaw = [PreProcessedData,filesep,Prefix,filesep,'tr2dProject',filesep,'RAW.tif'];
if exist( filenameRaw )
    Answer=input('Nuclear data has already been exported for tr2d. Do you want to export again? (y/N)','s');
    if strcmpi(Answer,'y')
        delete(filenameRaw);
    else
        CreateStack=0;
    end
end

if CreateStack
    %Load all nuclear images
    D=dir([PreProcessedData,filesep,Prefix,filesep,'*-His_*.tif']);
    for i=1:length(D)
        Image=imread([PreProcessedData,filesep,Prefix,filesep,...
            D(i).name]);
        imwrite(Image,filenameRaw,'WriteMode','append')
    end
end

%Pull out the pixel size (in um) and the list of time stamps (in seconds).
%Save it as a CSV file that will serve a input for tr2d.
load([DynamicsResults,filesep,Prefix,filesep,'FrameInfo.mat'])
PixelSize=FrameInfo(1).PixelSize;
TimeStamps=[FrameInfo.Time];
%Save the information as CSV
CSVOutput(2,:)=TimeStamps;
SCSVOutput(1,1)=PixelSize;
csvwrite([PreProcessedData,filesep,Prefix,filesep,'tr2dProject',...
    filesep,'mRNADynamics.csv'],CSVOutput);


%%%%%%%
%RUN TR2D
%%%%%%%%
% Check if IJM is started already and start if it is not
% try
%     %this is just some function that can only be called if IJM is set up
%     IJM.getIdentifier() 
% catch
%     addpath('e:/Fiji.app/scripts') % Update for your ImageJ installation
%     ImageJ                         % Initialize IJM and MIJ
% end
% Start tr2d
%MIJ.run('Tr2d 0.2.2-SNAP');
%params = javaArray('java.lang.String', 1);
%params(1) = java.lang.String(['-p ',PreProcessedData,filesep,Prefix,filesep,'tr2dProject'])
%tr2d = com.indago.tr2d.app.garcia.Tr2dApplication()
%tr2d.main(params)
%MIJ.exit()
scriptPath = fileparts(mfilename('fullpath'));
jar = fullfile(scriptPath, 'tr2d.jar');
projectPath = [PreProcessedData,filesep,Prefix,filesep,'tr2dProject'];
exportPath = [projectPath,filesep,'mRNADynamicsExport'];
system(['java -jar ', jar, ' -p ', projectPath, ' -e ', exportPath]);

%Check that we have the tr2d results and import them
if exist([PreProcessedData,filesep,Prefix,filesep,'tr2dProject',filesep,...
    'mRNADynamicsExport',filesep,'tr2d_objects.csv'], 'file')&&...
    exist([PreProcessedData,filesep,Prefix,filesep,'tr2dProject',filesep,...
    'mRNADynamicsExport',filesep,'tr2d_tracks.csv'], 'file')

    ImportTr2d(Prefix)
else
    error('No tr2d results found, did you export?')    
end








