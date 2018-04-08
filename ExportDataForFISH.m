function Prefix=ExportDataForFISH(varargin)

%This function grabs individual z-stacks and splits them in
%multiple channels so that it can be analyzed by the FISH code.
%It adds a blank image at the beginning and the end so that the code
%doesn't discard columns that peak at the edges of the Z-stack. 

%Options:
%medianprojection: Uses a median projection in the nuclear channel rather
%                  than the default maximum projection


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


%Parameters:
NIndices=3;     %Number of indices ScanImage used to save the files
MaxShift=9;     %Maximum shift in pixels corresponding to image shift and
                %alignment
MaxHistone=1000;    %Maximum intensity for the histone channel. Anything above
                    %this will be capped.
ProjectionType = 'maxprojection'; %Default setting for z-projection is maximum-based...
                    %This may fail when high intensity reflections are present
                

%Look at parameters
PrefixOverrideFlag = 0;
SkipFrames=[];
k=1;
while k<=length(varargin)
    if strcmpi(varargin{k},'skipframes')
        SkipFrames=varargin{k+1};
        k=k+1;
        warning('SkipFrame mode.')
    elseif strcmpi(varargin{k},'medianprojection')
        ProjectionType = 'medianprojection';
    else
        Prefix = varargin{k};
        PrefixOverrideFlag = 1;
    end
    k=k+1;
end

[SourcePath,FISHPath,DropboxFolder,MS2CodePath, PreProcPath,...
    Folder, Prefix, ExperimentType, Channel1, Channel2,OutputFolder, Channel3...
    ] = readMovieDatabase(PrefixOverrideFlag);

[D, FileMode] = DetermineFileMode(Folder);

%Create the output folder
OutputFolder=[PreProcPath,filesep,Prefix];
mkdir(OutputFolder)

%Generate FrameInfo
FrameInfo=struct('LinesPerFrame',{},'PixelsPerLine',{},...
    'NumberSlices',{},'ZStep',{},'FileMode',{},...
    'PixelSize',{});


%Extract channel information
%This information will be stored in FrameInfo for use by subsequent parts
%of the code. Note, however, that the channels are also extracted in this
%code for each data type. I should integrate this.




if strcmp(FileMode,'TIF') && ~strcmp(FileMode,'DSPIN')

  FrameInfo = process2PhotonPrincetonData(Folder, D, FrameInfo, Channel2, OutputFolder);

elseif strcmp(FileMode, 'LAT')
      
  FrameInfo = processLatticeLightSheetData(Folder, D, Channel1, Channel2, ProjectionType, Prefix, OutputFolder);
                  
%LSM mode
elseif strcmp(FileMode,'LSM')
    
  FrameInfo = processZeissConfocalLSMData(Folder, D, ExperimentType, Channel1, Channel2, Prefix, OutputFolder);
    
%LIFExport mode
elseif strcmp(FileMode,'LIFExport')
    
  FrameInfo = processLIFExportMode(Folder, ExperimentType, ProjectionType, Channel1, Channel2, Channel3, Prefix, OutputFolder);        
    
%Nikon spinning disk confocal mode - TH/CS 2017
elseif strcmp(FileMode,'DSPIN')||strcmp(FileMode,'DND2')

  FrameInfo = processSPINandND2Data(Folder, D, ExperimentType, Channel1, Channel2, SourcePath, Prefix, OutputFolder, DropboxFolder);
  
end

doFrameSkipping(SkipFrames, FrameInfo, OutputFolder);

%Save the information about the various frames
mkdir([DropboxFolder,filesep,Prefix])
save([DropboxFolder,filesep,Prefix,filesep,'FrameInfo.mat'],...
    'FrameInfo')
