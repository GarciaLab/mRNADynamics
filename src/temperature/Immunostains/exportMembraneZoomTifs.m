function exportMembraneZoomTifs(Prefix, varargin)
%%
%
%   function exportImmunostainTIFStacks(AllImages, imagingModality, NChannels, NReplicates, NSlices,...
%   Prefix, moviePrecision, hisPrecision)





[rawDataPath,~,~,~, ~, ~, movieDatabasePath, movieDatabaseFolder]=...
    DetermineLocalFolders;
Subfolder = [Prefix(1:10),filesep,Prefix(12:length(Prefix))];
rawDataFolder = strcat(rawDataPath,filesep,Subfolder);
MembranePath = [rawDataFolder, filesep, 'Membrane', filesep, 'MembraneZoomInfo.lif'];
MembraneFolder = [rawDataFolder, filesep, 'Membrane'];

if isfile(MembranePath)
    
    liveExperiment = LiveExperiment(Prefix);
    PreProcFolder = liveExperiment.preFolder;
    
    [XMLFolder, seriesPropertiesXML, seriesXML] = getSeriesFiles(MembraneFolder);
    
    [LIFImages, LIFMeta] = loadLIFFile(MembraneFolder);
    if ~isempty(str2double(LIFMeta.getPixelsPhysicalSizeX(0))) &...
            ~isnan(str2double(LIFMeta.getPixelsPhysicalSizeX(0)))
        PixelSize_um = str2double(LIFMeta.getPixelsPhysicalSizeX(0));
    else
        try
            PixelSize_um = str2double(LIFMeta.getPixelsPhysicalSizeX(0).value);
        catch %no idea man
            PixelSize_um = str2double(LIFMeta.getPixelsPhysicalSizeX(1).value);
        end
    end
    
    ySize = size(LIFImages{1}{1,1}, 1);
    xSize = size(LIFImages{1}{1,1}, 2);
    zSize = size(LIFImages{1}, 1);
    NEmbryos = size(LIFImages, 1);
    
    memMat = zeros(ySize, xSize,zSize, NEmbryos, 'uint16');
    
    NameSuffix = '_ZoomMembrane';
    % loop through each series
    for embryoIndex = 1:NEmbryos
        NewName = [Prefix, '_Position',iIndex(embryoIndex,3),...
                NameSuffix, '.tif'];
            
        for zIndex = 1:zSize
            
            
            
            memMat(:,:,zIndex, embryoIndex) = LIFImages{embryoIndex,1}{zIndex,1};
            if zIndex == 1
                imwrite(squeeze(memMat(:,:,zIndex, embryoIndex)), [PreProcFolder, filesep, NewName]);
            else
                imwrite(squeeze(memMat(:,:,zIndex, embryoIndex)), [PreProcFolder, filesep, NewName],...
                     'WriteMode', 'append');
            end
        end
        
        
        
    end

    MemZoomOutPath = [PreProcFolder, filesep, NewName];
    
    MembraneZoomResultsFolder = [liveExperiment.resultsFolder, filesep, 'ZoomMembraneInfo'];
    if ~exist(MembraneZoomResultsFolder, 'dir')
        mkdir(MembraneZoomResultsFolder);
        
    end
    
    MembraneZoomPixelPath = [MembraneZoomResultsFolder, filesep, 'MembraneZoomPixelSize.mat'];
    save(MembraneZoomPixelPath,'PixelSize_um');
    
    MembraneZoomImageSizePath = [MembraneZoomResultsFolder, filesep, 'MembraneZoomImageSize.mat'];
    save(MembraneZoomImageSizePath,'xSize', 'ySize', 'zSize');
   
    medMemMat = zeros(ySize, xSize, NEmbryos, 'uint16'); 
    for embryoIndex = 1:NEmbryos
 
            
            
            
            medMemMat(:,:, embryoIndex) = median(squeeze(memMat(:,:,:,embryoIndex)),3);
        
        
        
        
    end
    MedNameSuffix = 'MedZoomMembrane';
    MedNewName = [Prefix,'-',MedNameSuffix, '.tif'];
    
    
    MedMemZoomOutPath = [PreProcFolder, filesep, MedNewName];
    saveNuclearProjection(medMemMat, MedMemZoomOutPath);
end






