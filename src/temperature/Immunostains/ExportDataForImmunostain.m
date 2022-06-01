% ExportDataForImmunostain([Options])
%
% DESCRIPTION
% %This function grabs individual z-stacks and splits them in
% multiple channels so that it can be analyzed by the FISH code.
% It adds a blank image at the beginning and the end so that the code
% doesn't discard columns that peak at the edges of the Z-stack.
%
% ARGUMENTS
%
% [Options]: See below.
%
% OPTIONS
%
% Prefix: If you want to avoid the user input selection of the data folder,
%           you can pass your desired Prefix as the first option. It MUST 
%           be the first option, and you CANNOT use the 'rootFolder'
%           option with this option, otherwise it will break.
%           Prefix standard format is 'YYYY-MM-DD-YOURPROJECTNAME'.
%           E.g. ExportDataForLivemRNA('2020-02-03-Dl_Ven_snaBAC_mCh_01',
%                                       'medianprojection')
%
% 'medianprojection': Uses a median projection in the nuclear channel rather
%                  than the default maximum projection
% 'middleprojection': Uses a max projection in the nuclear channel, but only
%                   of the middle slices (11-16) to prevent bright
%                   reflections from overpowering the signal
% 'meanprojection'  : Uses the mean projection in the nuclear channel
%                     instead of the default maximum projection. 
% 'PreferredFileForTest': Uses the given filename for the prefered FF file.
%                         PreferredFileForTest needs to be a struct with 
%                         a field name called "fileName" with the file name. 
%                         EL 12/3/2018
% 'skipframes',framesToSkip: deletes frames given by framesToSkip
%                            EL 12/3/2018 option might not be working as
%                            intended. 
% 'keepTifs': MPF 9/12/2018 Do not delete source folder TIF files when
% running. This is used for testing purposes.
% 'generateTifs': MPF 11/11/2018 Additionally run filterMovie to generate Tifs stacks
% 'skipExtraction': This doesn't extract LIF files to Tifs. Occasionally
%                   useful if only the FrameInfo is desired. 
% 'rootFolder',rootFolder: open a directory different from the default user directory
%               if there is no given prefix. This is useful if you are
%               opening data that is in a different project folder.
% 'zslicesPadding': if series have different number of z-slices, pad them
% with blank images so every generates series has the same amount
% 'nuclearGUI': accepts true (default) or false if you want to open the
% nuclear channel / anaphase frame choosing GUI 
% 'skipNuclearProjection': only extract channels 
% 
% OUTPUT
% Exported tif images are placed in the PreProcessedData folder and divided
% into a nuclear channel for tracking/protein quantification and a channel
% for transcriptional loci.
% All TIFs from source folder are deleted except argument 'keepTifs' is
% used.
%
% Author (contact): Hernan Garcia
% Created: 01/01/2016
% Last Updated: 10/13/2019
%
% Documented by: Armando Reimer (areimer@berkeley.edu)
%
%Note:
%Flatfield image: We assume there is no background of fluorescence / dark
%current.
%
%
%Where do things go?
%
%1) The original raw data: RawDynamics Folder
%2) The data processed such that it can be analyzed by the FISH code:
%   the PreProcessedData Folder
%3) The data analyzed by the FISH code: ProcessedData folder
%4) The resulting structures from the particle tracking:
%   DynamicsResults folder
%
%The idea of (4) being in Dropbox is that I don't need to be synchronizing
%the part related to the manual analysis.
% 
function Prefix = ExportDataForImmunostain(varargin)

cleanupObj = onCleanup(@myCleanupFun);
clear LiveExperiment;

warning('off', 'MATLAB:MKDIR:DirectoryExists');





  [Prefix, SkipFrames, ProjectionType, PreferredFileNameForTest, ~,...
    generateTifStacks, nuclearGUI, skipExtraction, rootFolder, zslicesPadding,...
    dataType, skipNuclearProjection] = ...
    ...
    exportDataForLivemRNA_processInputParameters(varargin{:});

  [rawDataPath, ProcPath, DropboxFolder, ~, PreProcPath, rawDataFolder, Prefix, ExperimentType, Channel1, Channel2, ~,...
    Channel3, ~, ~, ~, Channel4, Channel5] = readMovieDatabase(Prefix,'rootFolder', rootFolder);

Channels = {Channel1, Channel2, Channel3, Channel4, Channel5};
%%

%   if ~isempty(rootFolder)
    [D, FileMode] = DetermineFileMode(rawDataFolder);
%   else
%     [D, FileMode] = DetermineFileMode(rootFolder);
%   end

%Create the output folder
OutputFolder = [PreProcPath, filesep, Prefix];
disp(['Creating folder: ', OutputFolder]);
mkdir(OutputFolder)
  
%Create the results folder
DropboxFolderName = [DropboxFolder, filesep, Prefix];
disp(['Creating folder: ', DropboxFolderName]);
mkdir(DropboxFolderName);

load('ReferenceHist.mat', 'ReferenceHist');

  %Generate FrameInfo
  FrameInfo = struct('LinesPerFrame', {}, 'PixelsPerLine', {}, ...
    'NumberSlices', {}, 'ZStep', {}, 'FileMode', {}, ...
    'PixelSize', {});

 %%
 liveExperiment = LiveExperiment(Prefix);
disp('Exporting movie file...');

resultsFolder = liveExperiment.resultsFolder;

moviePrecision = 'uint16';
hisPrecision = 'uint16';


%Load the reference histogram for the fake histone channel
% load('ReferenceHist.mat', 'ReferenceHist');

markandfind = true;

[XMLFolder, seriesPropertiesXML, seriesXML] = getSeriesFiles(rawDataFolder);


try
    if contains(seriesPropertiesXML(1).name, 'Mark_and_Find')
        markandfind = true;
    end
catch % do nothing
end

[LIFImages, LIFMeta] = loadLIFFile(rawDataFolder);
 %Obtains frames information
 [NSeries, NFrames, NSlices,...
     NPlanes, NChannels, Frame_Times] = getFrames(LIFMeta, ExperimentType);
 
 %use the old method(exported from lasx) if the files are exported
 %already. if they're not, just use bioformats. the lasx method is being
 %deprecated.
 if ~isempty(XMLFolder)
     timeStampRetrievalMethod = 'lasx';
 else
     timeStampRetrievalMethod = 'bioformats';
 end
 
 if sum(NFrames)~=0
     
     switch timeStampRetrievalMethod
         
         
         case 'manual'
             
             xml_file = [liveExperiment.rawFolder, filesep, 'lifMeta.xml'];
             
             if ~exist(xml_file, 'file')
                 generateLIFMetaDataXML(Prefix, xml_file);
             end
             
             Frame_Times = getTimeStampsFromLifXML(xml_file);
             
         case 'lasx'
             
             Frame_Times = obtainFrameTimes(XMLFolder, seriesPropertiesXML,...
                 NSeries, NFrames, NSlices, NChannels);
             
         case 'bioformats'
             
             Frame_Times = getFrameTimesFromBioFormats(LIFMeta, NSlices);
             
         otherwise, error('what?');
             
     end
 end
 
 [InitialStackTime, zPosition] = getFirstSliceTimestamp(NSlices,...
     NSeries, NPlanes, NChannels, Frame_Times, XMLFolder, seriesXML);
 
 FrameInfo = recordFrameInfo(NFrames, NSlices, InitialStackTime, LIFMeta, zPosition);
 
 if markandfind
     FrameInfo = repmat(FrameInfo, NSeries, 1);
 end
 
 save([resultsFolder, filesep, 'FrameInfo.mat'], 'FrameInfo', '-v6');
    

 if sum(NFrames) == 0
     NFrames = ~NFrames;
 end
    
%Obtains frames information
[NSeries, NSlices,NPlanes,...
    NChannels, NReplicates] = getMarkAndFindInfo(LIFMeta, ExperimentType);

MarkAndFindInfo.NSeries = NSeries;
MarkAndFindInfo.NSlices = NSlices;
MarkAndFindInfo.NPlanes = NPlanes;
MarkAndFindInfo.NChannels = NChannels;
MarkAndFindInfo.NReplicates = NReplicates;
  save([DropboxFolder, filesep, Prefix, filesep, 'MarkAndFindInfo.mat'], 'MarkAndFindInfo', '-v6');

%use the old method(exported from lasx) if the files are exported
%already. if they're not, just use bioformats. the lasx method is being
%deprecated.
%Find the flat field (FF) information
LIFExportMode_flatFieldImage(LIFMeta,...
    rawDataFolder, PreProcPath, Prefix, PreferredFileNameForTest);



% this function exports tif z stacks

exportImmunostainTIFStacks(LIFImages, 'LIF', NChannels, NReplicates, NSlices,...
    Prefix, moviePrecision, hisPrecision);
% try
%     exportMembraneZoomTifs(Prefix);
% catch
%     disp('Missing Membrane Zoom Tifs')
% end

chooseIHProjections(...
    Prefix,'ProjectionType',ProjectionType,...
    'ReferenceHist',ReferenceHist);


    



disp('Movie files exported.');

 
 
 
 %%
 
 %make folders we'll need later
 DogOutputFolder=[ProcPath,filesep,Prefix, '_', filesep, 'dogs',filesep];
 mkdir(DogOutputFolder);
    

%end
