function [FileMode, EmbryoName, projectDate] = getMicroscope(Prefix)

RawDynamicsPath= DetermineAllLocalFolders(Prefix)

%Determine whether we're dealing with 2-photon data from Princeton or LSM
%data. 2-photon data uses TIF files. In LSM mode multiple files will be
%combined into one.
%Find out the date it was taken
Dashes=strfind(Prefix,'-');
projectDate=Prefix(1:Dashes(3)-1);
EmbryoName=Prefix(Dashes(3)+1:end);
rawPrefixPath = [RawDynamicsPath,filesep,projectDate,filesep,EmbryoName,filesep];


DTIF=[dir([rawPrefixPath,'*.tif']), dir([rawPrefixPath,'*.jpg'])];
DLSM=dir([rawPrefixPath,'*.lsm']);
DCZI=dir([rawPrefixPath,'*.czi']);
DLIF=dir([rawPrefixPath,'*.lif']);
DLAT=dir([rawPrefixPath,'*_Settings.txt']);
DSPIN=dir([rawPrefixPath,'FullEmbryo',filesep,'*.nd']);     %Nikon spinning disk . CS20170911

% OME-TIFF xml companion file (*.ome).
% Not mandatory per OME-TIFF standard, but our examples have it so for know we detect ome-tiff based on the presence of this file.
OMETIFF = dir([rawPrefixPath,'*.ome']);

if ~isempty(OMETIFF)
    disp('OME-TIFF with .ome XML companion file mode')
    FileMode = 'OMETIFF';
elseif ~isempty(DTIF) & isempty(DLSM)& isempty(DSPIN)
    if isempty(DLIF)
        if isempty(DLAT)
            disp('2-photon @ Princeton data mode')
            FileMode='TIF';
        else
            disp('Lattice Light Sheet data mode')
            FileMode='LAT';
        end
    else
        disp('LIF export mode')
        FileMode='LIFExport';
    end
elseif (isempty(DTIF))&(~isempty(DLSM))
    disp('LSM mode')
    FileMode='LSM';
elseif (isempty(DTIF))&(~isempty(DCZI))
    disp('CZI mode')
    FileMode='CZI';
elseif (~isempty(DSPIN))        %CS20170911
    disp('Nikon spinning disk mode with .nd files')
    FileMode='DSPIN';
else
    error('File type not recognized')
end

end