% ExportDataForLivemRNA([Options])
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
% 'medianprojection': Uses a median projection in the nuclear channel rather
%                  than the default maximum projection
% 'middleprojection': Uses a max projection in the nuclear channel, but only
%                   of the middle slices (11-16) to prevent bright
%                   reflections from overpowering the signal
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
% Last Updated: 8/30/2018
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

function Prefix = ExportDataForLivemRNA(varargin)

  [Prefix, SkipFrames, ProjectionType, PreferredFileNameForTest, keepTifs,...
    generateTifs] = exportDataForLivemRNA_processInputParameters(varargin{:});

  [SourcePath, ~, DropboxFolder, ~, PreProcPath, Folder, Prefix, ExperimentType, Channel1, Channel2, ~,...
    Channel3] = readMovieDatabase(Prefix);

  [D, FileMode] = DetermineFileMode(Folder);

  %Create the output folder
  OutputFolder = [PreProcPath, filesep, Prefix];
  mkdir(OutputFolder)

  %Generate FrameInfo
  FrameInfo = struct('LinesPerFrame', {}, 'PixelsPerLine', {}, ...
    'NumberSlices', {}, 'ZStep', {}, 'FileMode', {}, ...
    'PixelSize', {});

  %Extract channel information
  %This information will be stored in FrameInfo for use by subsequent parts
  %of the code. Note, however, that the channels are also extracted in this
  %code for each data type. I should integrate this.
  if strcmpi(FileMode, 'TIF')
    %Maximum shift in pixels corresponding to image shift and alignment
    MaxShift = 9;

    %Maximum intensity for the histone channel. Anything above this will be capped.
    MaxHistone = 1000;

    FrameInfo = process2PhotonPrincetonData(Folder, D, FrameInfo, Channel2, MaxShift, MaxHistone, OutputFolder);
  elseif strcmpi(FileMode, 'LAT')
    FrameInfo = processLatticeLightSheetData(Folder, D, Channel1, Channel2, ProjectionType, Prefix, OutputFolder);

  elseif strcmpi(FileMode, 'LSM')
    FrameInfo = processZeissConfocalLSMData(Folder, D, FrameInfo, ExperimentType, Channel1, Channel2, ProjectionType, Prefix, OutputFolder);

  elseif strcmpi(FileMode, 'LIFExport')
    FrameInfo = processLIFExportMode(Folder, ExperimentType, ProjectionType, Channel1, Channel2, Channel3, Prefix, ...
      OutputFolder, PreferredFileNameForTest, keepTifs);

  elseif strcmpi(FileMode, 'DSPIN') || strcmpi(FileMode, 'DND2')
    %Nikon spinning disk confocal mode - TH/CS 2017
    FrameInfo = processSPINandND2Data(Folder, D, FrameInfo, ExperimentType, Channel1, Channel2, SourcePath, Prefix, OutputFolder, DropboxFolder);
  end

  doFrameSkipping(SkipFrames, FrameInfo, OutputFolder);

  %Save the information about the various frames
  mkdir([DropboxFolder, filesep, Prefix]);
  save([DropboxFolder, filesep, Prefix, filesep, 'FrameInfo.mat'], 'FrameInfo');

  if generateTifs
    filterMovie(Prefix, 'Tifs');
    disp(['Prefix: ', Prefix]);
  end

end