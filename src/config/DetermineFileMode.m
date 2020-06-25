%Determine whether we're dealing with 2-photon data from Princeton or LSM
%data. 2-photon data uses TIF files. In LSM mode multiple files will be
%combined into one.
function [D, FileMode] = DetermineFileMode(rawDataFolder)

DTIF=dir([rawDataFolder,filesep,'*.tif']);
DLSM=dir([rawDataFolder,filesep,'*.lsm']);     %Zeiss confocal, old LSM format
DLIF=dir([rawDataFolder,filesep,'*.lif']);     %Leica confocal
DCZI=dir([rawDataFolder,filesep,'*.czi']);     %Zeiss confocal, new CZI format
DLAT=dir([rawDataFolder,filesep,'*_Settings.txt']);
DSPIN=dir([rawDataFolder,filesep,'*.nd']);     %Nikon spinning disk
DND2=dir([rawDataFolder,filesep,'*.nd2']);    %Nikon point scanner .nd2 files
OMETIFF = dir([rawDataFolder,filesep,'*.ome']);
% OME-TIFF xml companion file (*.ome).
% Not mandatory per OME-TIFF standard, but our examples have it so for know we detect ome-tiff based on the presence of this file.

if ~isempty(OMETIFF)
    disp('OME-TIFF with .ome XML companion file mode')
    D = OMETIFF;
    FileMode = 'OMETIFF';
elseif ~isempty(DLIF)
    disp('LIF export mode')
    D='';
    FileMode='LIFExport';
elseif ~isempty(DTIF) && isempty(DLSM)...
        && isempty(DCZI) && isempty(DSPIN) && isempty(DLIF)
    
    if isempty(DLAT)
        disp('2-photon @ Princeton data mode')
        D=DTIF;
        FileMode='TIF';
    else
        disp('Lattice Light Sheet data mode')
        D=DTIF;
        FileMode='LAT';
    end
elseif isempty(DTIF) && ~isempty(DLSM)
    disp('LSM mode')
    D=DLSM;
    FileMode='LSM';
elseif isempty(DTIF) && ~isempty(DCZI)
    disp('LSM (CZI) mode')
    D=DCZI;
    FileMode='LSM';
elseif ~isempty(DSPIN)
    disp('Nikon spinning disk mode with .nd files')
    if isempty(DTIF)
        error('No TIF files found')
    end
    D=DTIF;
    FileMode = 'DSPIN';
elseif ~isempty(DND2)
    disp('Nikon LSM Mode');
    D=DND2;
    FileMode='DND2';
else
    error('File type not recognized.')
end


end