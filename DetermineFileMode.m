%Determine whether we're dealing with 2-photon data from Princeton or LSM
%data. 2-photon data uses TIF files. In LSM mode multiple files will be
%combined into one.
function [D, FileMode] = DetermineFileMode(Folder)
    DTIF=dir([Folder,filesep,'*.tif']);
    DLSM=dir([Folder,filesep,'*.lsm']);     %Zeiss confocal, old LSM format
    DLIF=dir([Folder,filesep,'*.lif']);     %Leica confocal
    DCZI=dir([Folder,filesep,'*.czi']);     %Zeiss confocal, new CZI format
    DLAT=dir([Folder,filesep,'*_Settings.txt']);
    DSPIN=dir([Folder,filesep,'*.nd']);     %Nikon spinning disk
    DND2=dir([Folder,filesep,'*.nd2']);    %Nikon point scanner .nd2 files

    if ~isempty(DTIF) && isempty(DLSM) && isempty(DCZI) && isempty(DSPIN)
        if isempty(DLIF)
            if isempty(DLAT)
                disp('2-photon @ Princeton data mode')
                D=DTIF;
                FileMode='TIF';
            else
                disp('Lattice Light Sheet data mode')
                D=DTIF;
                FileMode='LAT';
            end
        else
            disp('LIF export mode')
            D=DTIF;
            FileMode='LIFExport';
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
        error('File type not recognized. For LIF files, were they exported to TIF?')
    end
end