function [CytoFluoDensity] = processLlamaTagData_subtractBGFluo(varargin)
% This function calculates the average fluorescence intensity of the whole cytoplasm for
% each frame and z slice. 
% It does so by first creating a binary mask using the nuclear segmentation information stored in Ellipses.mat
% it then applies this mask to the PreProcessed Nuclear Protein (:input) files

% Input: A movie prefix

% Output: a [frame x Z] array containing the average pixel intensity of the
% cytoplasm for each z slice in each frame.

% Authors: Simon Alamos and Jordan Xiao (modified by YangJoon Kim)
% Contact: simon.alamos@berkeley.edu

% November 2017
%% Get Folder info, set up stuff, etc

[SourcePath,FISHPath,DefaultDropboxFolder,MS2CodePath,PreProcPath]=...
    DetermineLocalFolders;

Prefix=varargin{1};

FilePrefix=[Prefix,'_'];

%% Now get the actual Dropbox folder
[SourcePath,FISHPath,DropboxFolder,MS2CodePath,PreProcPath]=...
    DetermineLocalFolders(Prefix);

% This is to get the name of the channels from MovieDatabase
[XLSNum,XLSTxt,XLSRaw]=xlsread([DropboxFolder,filesep,'MovieDatabase.csv']);
Channel1Column=find(strcmp(XLSTxt(1,:),'Channel1'));
Channel2Column=find(strcmp(XLSTxt(1,:),'Channel2'));
DataFolderColumn=find(strcmp(XLSRaw(1,:),'DataFolder'));
Dashes=findstr(Prefix,'-');
PrefixRow=find(strcmp(XLSRaw(:,DataFolderColumn),[Prefix(1:Dashes(3)-1),'\',Prefix(Dashes(3)+1:end)]));
if isempty(PrefixRow)
    PrefixRow=find(strcmp(XLSRaw(:,DataFolderColumn),[Prefix(1:Dashes(3)-1),'/',Prefix(Dashes(3)+1:end)]));
    if isempty(PrefixRow)
        error('Could not find data set in MovieDatabase.XLSX. Check if it is defined there.')
    end
end

Channel1=XLSTxt(PrefixRow,Channel1Column);
Channel2=XLSTxt(PrefixRow,Channel2Column);

%Load all the information

%load([DropboxFolder,filesep,Prefix,filesep,'Particles.mat'])
%load([DropboxFolder,filesep,Prefix,filesep,'CompiledParticles.mat'])
%load([DropboxFolder,filesep,Prefix,filesep,'CompiledNuclei.mat'])
load([DropboxFolder,filesep,Prefix,filesep,[Prefix '_lin.mat']])
load([DropboxFolder,filesep,Prefix,filesep,'Ellipses.mat'])


%Check that FrameInfo exists
if exist([DropboxFolder,filesep,Prefix,filesep,'FrameInfo.mat'])
    load([DropboxFolder,filesep,Prefix,filesep,'FrameInfo.mat'])
end

% %get frames of each mitosis
[XLSNum,XLSTxt,XLSRaw]=xlsread([DefaultDropboxFolder,filesep,'MovieDatabase.xlsx']);
% nc10Column=find(strcmp(XLSRaw(1,:),'nc10'));
% nc11Column=find(strcmp(XLSRaw(1,:),'nc11'));
% nc12Column= strcmp(XLSRaw(1,:),'nc12');
% nc13Column=find(strcmp(XLSRaw(1,:),'nc13'));
% nc14Column=find(strcmp(XLSRaw(1,:),'nc14'));

DataFolderColumn=find(strcmp(XLSRaw(1,:),'DataFolder'));
Dashes=findstr(Prefix,'-');
PrefixRow=find(strcmp(XLSRaw(:,DataFolderColumn),[Prefix(1:Dashes(3)-1),'\',Prefix(Dashes(3)+1:end)]));

%% Create a 3D mask using ellipses

% Our nucleus segmentation algorithm  DOES NOT track nuclei in Z. In
% consequence, a single segmentation is applied to all Z slices in any
% given frame. This means we don't need 4 dimensions, but just 3: X,Y,T

% create a black image of the same X,Y size as the data

%get dimensions
TotalFrames = length(FrameInfo);
Rows = FrameInfo(1).LinesPerFrame;
Columns = FrameInfo(1).PixelsPerLine;

%create black image
BlackImage = zeros(Rows,Columns);

%create empty XYT array to store masks
MovieMask = ones(Rows,Columns,TotalFrames);

% Define the scaling factor for the nuclear mask
ScaleFactor = 1.4 % 40% increase

figure(1)
h = waitbar(0,'Please wait, creating a XYT 3D mask for the cytoplasm');
for frame = 1:TotalFrames
    
    Mask = BlackImage; %we will add this frame's ellipses to this black image
    
    FrameEllipses = Ellipses{frame}; %get the ellipses belonging to this frame
    CentersX = FrameEllipses(:,1); %get the centers of this frame's ellipses
    CentersY = FrameEllipses(:,2); %get the sizes of this frame's ellipses
    EllipsesRadius = FrameEllipses(1,3); %get the radius of the ellipses in this frame
    
    %now we will add the ellipses one by one as white circles
    for ellipse = 1:size(FrameEllipses,1)
        
        % get parameters needed to create a disk
        SingleEllipseMask = zeros(Rows,Columns);
        EllipseXCenter = CentersX(ellipse);
        EllipseYCenter = CentersY(ellipse);
        % Radius of the ellipse. 
        % Increase with scaling factor to make sure we get the whole nucleus
        Radius = ScaleFactor*(ceil(EllipsesRadius));
        
        
        % create a white disk where the ellipse is
        Circle=BlackImage;
        Circle=MidpointCircle(Circle,Radius,EllipseYCenter,...
        EllipseXCenter,1); 
    
        % add disk (ellipse) to the mask
        Mask = Mask+Circle; 

    end
    
    Mask = Mask==0; % invert the mask
    MovieMask(:,:,frame) = Mask;  %add mask to the XYT movie array of masks
    imshow(Mask,[])
    %pause(0.01)
    waitbar(frame / TotalFrames);


end
close(h)
%% Apply mask to the Protein channel (movie)

% First, we will create a 4D movie using the Bcd / Dorsal channel so that we can
% apply the mask to it by doing element-wise matrix multiplication.

%get the size of the Z dimension
Slices = FrameInfo(1).NumberSlices + 2;

% Figure out which channel has Protein (:input)
if contains(lower(Channel1{1}),'input') || contains(lower(Channel1{1}),'bcd')||contains(lower(Channel1{1}),'dl')
    InputChannel = 'ch01';
elseif contains(lower(Channel2{1}),'input') || contains(lower(Channel2{1}),'bcd')||contains(lower(Channel2{1}),'dl')
    InputChannel = 'ch02';
else
    display ('no channel with Bcd or Dl in its name found in this dataset')
end

%InputChannel = 'ch01';

InputChannelFiles = dir([PreProcPath,filesep,Prefix,filesep,'*' InputChannel '.tif']); % save info about files in a struct
Frames = length(InputChannelFiles)/Slices; % get the number of frames

%%  We need add to the struct the Z and T info of each file, which is stored in their names

for i= 1:length(InputChannelFiles)
    
    ImageName = InputChannelFiles(i).name; %get the name string
    SplitName = strsplit(ImageName,'_'); %split the name in strings between '_' characters
    
    % Z information
    ZSlice1 = SplitName{end-1}; %select only the part of the name referring to the Z info
    ZSlice2 = str2num(ZSlice1(2:end)); %get an actual number corresponding to the z slice
 
    % Time information
    Time1 = SplitName{end-2}; %select only the part of the name referring to the T info
    Time2 = str2num(Time1); %get an actual number corresponding to the frame
    
    % Add info to the struct
    InputChannelFiles(i).Z = ZSlice2;
    InputChannelFiles(i).Time = Time2;
    
end

%% Now, loop over time and z to apply the cytoplasm mask to each image
% We will store the cytoplasm fluorescence in a T x Z 2D matrix called
% 'CytoFluo'
% For each frame we will calculate the median cytoplasmid fluorescence over
% Z-stacks, and also averaged over binned area.

CytoFluoImage = zeros(Rows,Columns,Slices);
CytoFluoDensity = zeros(Rows,Columns,Frames); %initialize the array where we'll store the data
MovieMask(MovieMask==0)=NaN; %convert 0s to NaNs because the raw data contains 0s

D = [PreProcPath,filesep,Prefix];
for T = 1:Frames
    Mask = MovieMask(:,:,T);
    
    for Z = 1:Slices
        
        for i = 1:length(InputChannelFiles)
            
            ImageName = InputChannelFiles(i).name;
            ZSlice = InputChannelFiles(i).Z;
            Time = InputChannelFiles(i).Time;
            
            if ZSlice == Z && Time == T
                
                NucleusFluoImage = imread([D,filesep,ImageName]); % read the image
                % imultiply sometimes doesn't work, depending on the
                % matlab version because of the class mismatch
                NucleusFluoImage = double(NucleusFluoImage);
                CytoFluoImage(:,:,Z) = immultiply(NucleusFluoImage,Mask); % apply mask to image
                %imshow(CytoFluoImage,[])
                %title([num2str(i) num2str(Z)])
                %waitforbuttonpress
                
            end
        end
    end
    % Calculate the median of cytofluo over Z-stacks
    % We can try maximum, but given the reflection in some of my movies, I
    % think this median makes more sense.
    CytoFluoDensity(:,:,T) = median(CytoFluoImage,3);
end


% save results in the DynamicsResults folder
save([DropboxFolder,filesep,Prefix,filesep,'CytoFluoDensity_median.mat'],'CytoFluoDensity')

%% Plot to check the masking (to see if the nuclei are included in the filter images)
% for i=1:Frames
%     imshow(CytoFluoDensity(:,:,i),[])
% %     pause
% end
% close all
%% Now let's actually caclulate the difference in fluorescence between nucleus and cytoplasm
%% (1) Calculate the averaged cyto fluo of whole field of view.
CytoFluoDensity_averaged = squeeze(nanmean(CytoFluoDensity,[1 2]));
CytoFluoDensity_std = squeeze(nanstd(CytoFluoDensity,0,[1 2]));
% Plot for checking
% errorbar(1:Frames, CytoFluoDensity_averaged, CytoFluoDensity_std)

% This is the cytoplasmic fluo density, per pixel, averaged for the whole
% area (field of view) at each frame. Since there's very little AP
% dependence (very little cytoplasmic proportion of TF), we will use this
% density to calculate the nuclear proportion of fluorophore.

% We need to get the integration area (in pixels) that was used to
% calculate the nuclear fluorescence density so that the units match with the cyto
% fluo denisty. This corresponds to the variable 'IntegrationRadius', which
% is defined in the function 'TrackNuclei' as being 2 (in microns)

absIntegrationRadius = 2; %in microns
pixelSize = FrameInfo(1).PixelSize; %in microns/pixel

% divide the standard two micrometers (assumed actual nucleus size) by the resolution ([pixelsize] = micrometers per pixel)
IntegrationRadius = floor(absIntegrationRadius/pixelSize); %in pixels
% not sure what this is for but it's in TrackNuclei so...
if ~mod(IntegrationRadius,2)
    IntegrationRadius=IntegrationRadius+1;
end

% This is taken from ExtractNuclearFluorescence. We use this to
% calculate the integration area.
Circle = logical(zeros(3*IntegrationRadius,3*IntegrationRadius));
Circle = MidpointCircle(Circle,IntegrationRadius,1.5*IntegrationRadius+0.5,...
    1.5*IntegrationRadius+0.5,1);
integrationArea = sum(sum(Circle)); %number of pixels in the integration area

% Multiply the integrationArea to CytoFluoDensity_averaged to get the
% cytoplasmic fluo proportion.
CytoFluo = CytoFluoDensity_averaged * integrationArea;

% Now, let's convert this cytoplasmic fluo to nuclear free FP fluo, by
% using the Supplementary from the LlamaTag paper.
% We will us the K_G of 0.8 from the paper.
Kg = 0.8; % Cyto to Nuclear free FP (eGFP)
NucFluo_freeFP = CytoFluo / Kg;

%% Loop over schnitzcells to calculate the delta (Nucleus-Cyto)
%for each schnitz, subtract the nuclear free FP fluorescence from the
%nuclear fluo, then assign a new field to those, such as
% schnitzcells.Fluo_BGsubtracted;
for schnitz = 1:length(schnitzcells)
    
    timeframes = schnitzcells(schnitz).frames; % vector with all the timeframes for a given schnitz
    schnitzFluo = schnitzcells(schnitz).Fluo; % array with one row per timeframe, one column per z slice
    
    NucFluoMax = max(schnitzFluo,[],2); % maximum fluo of all z-stacks (shoudl we try median?)

    % Subtract the Nuclear free FP from the NucFluoMax
    NucFluo = NucFluoMax - NucFluo_freeFP(timeframes);
    
    % Add the Nuc TF-eGFP fluo to the struct
    schnitzcells(schnitz).Fluo_BGsubtracted = NucFluo;
end

% Save this to schnitzcells
save([DropboxFolder,filesep,Prefix,filesep,Prefix,'_lin.mat'],'schnitzcells')
%% Plot to Check
% Check the Cytoplasmid fluo AP dependence.

% Plot the cyto fluo over X axis (AP)
% hold on
% for i=1:Frames
%     errorbar(1:Columns, movmean(nanmean(CytoFluoDensity(:,:,i),1),50)*integrationArea,nanstd(CytoFluoDensity(:,:,i),0,1)*integrationArea)
% end


%% (2) Binnig Cytoplasmic fluo using the AP information (in progress)
% % We're doing this so that we can calculate Averaged cytoplasmic fluo for
% % each AP bin which we can subtract from the nuclear fluo (MeanVectorAP).

%% (3) For each frame, Assign each non-zero pixel of cytofluo to the nearest schnitzcells.

end

    


