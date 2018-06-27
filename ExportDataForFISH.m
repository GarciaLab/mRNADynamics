%This function grabs individual z-stacks and splits them in
%multiple channels so that it can be analyzed by the FISH code.
%It adds a blank image at the beginning and the end so that the code
%doesn't discard columns that peak at the edges of the Z-stack. 

%Options:
%medianprojection: Uses a median projection in the nuclear channel rather
%                  than the default maximum projection
%middleprojection: Uses a max projection in the nuclear channel, but only
%                   of the middle slices (11-16) to prevent bright
%                   reflections from overpowering the signal


%Note:
%Flatfield image: We assume there is no background of fluorescence / dark
%current.


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
function Prefix = ExportDataForFISH(varargin)
addpath('LIFExport');
addpath('ZeissConfocalLSM');

[Prefix, SkipFrames, ProjectionType, PreferredFileNameForTest] = exportDataForFISH_processInputParameters(varargin{:})

[SourcePath,FISHPath,DropboxFolder,MS2CodePath, PreProcPath,...
    Folder, Prefix, ExperimentType, Channel1, Channel2,OutputFolder, Channel3...
    ] = readMovieDatabase(Prefix);

[D, FileMode] = DetermineFileMode(Folder);

%Create the output folder
OutputFolder = [PreProcPath,filesep,Prefix];
mkdir(OutputFolder)

%Generate FrameInfo
FrameInfo = struct('LinesPerFrame',{},'PixelsPerLine',{},...
    'NumberSlices',{},'ZStep',{},'FileMode',{},...
    'PixelSize',{});


%Extract channel information
%This information will be stored in FrameInfo for use by subsequent parts
%of the code. Note, however, that the channels are also extracted in this
%code for each data type. I should integrate this.
if strcmp(FileMode,'TIF')
  %Maximum shift in pixels corresponding to image shift and alignment
  MaxShift = 9; 

  %Maximum intensity for the histone channel. Anything above this will be capped.
  MaxHistone = 1000;

  FrameInfo = process2PhotonPrincetonData(Folder, D, FrameInfo, Channel2, MaxShift, MaxHistone, OutputFolder);
elseif strcmp(FileMode, 'LAT')
  FrameInfo = processLatticeLightSheetData(Folder, D, Channel1, Channel2, ProjectionType, Prefix, OutputFolder);

elseif strcmp(FileMode,'LSM')
  FrameInfo = processZeissConfocalLSMData(Folder, D, FrameInfo, ExperimentType, Channel1, Channel2, Prefix, OutputFolder);

elseif strcmp(FileMode,'LIFExport')
  FrameInfo = processLIFExportMode(Folder, ExperimentType, FrameInfo, ProjectionType, Channel1, Channel2, Channel3, Prefix, OutputFolder, PreferredFileNameForTest);        

elseif strcmp(FileMode,'DSPIN') || strcmp(FileMode,'DND2')
  %Nikon spinning disk confocal mode - TH/CS 2017
  FrameInfo = processSPINandND2Data(Folder, D, FrameInfo, ExperimentType, Channel1, Channel2, SourcePath, Prefix, OutputFolder, DropboxFolder);
end

doFrameSkipping(SkipFrames, FrameInfo, OutputFolder);

%Save the information about the various frames
mkdir([DropboxFolder, filesep, Prefix])
save([DropboxFolder, filesep, Prefix, filesep, 'FrameInfo.mat'], 'FrameInfo')
