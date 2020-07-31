function [MidImage,SurfImage,xShift,yShift] = rotateFullEmbryoSPIN(dirFullEmbryo, fullEmbryoPath, Channel1, Channel2)
%% 
% function [MidImage, SurfImage] = rotateFullEmbryoSPIN(dirFullEmbryo)
%
% DESCRIPTION
% Rotates the full embryo images and zoomed-in times series images to match
% each other for Nikon spinning disk data ('SPIN' file mode)
% (Might also do some other necessary full embryo image processing.)
%
% PARAMETERS
% dirFullEmbryo: directory for the folder containing the full embyro images
% fullEmbryoPath: path for the folder containing the full embyro images
% Channel1: first channel of the dataset
% Channel2: second channel of the dataset
%
% OPTIONS
% N/A
%
%
% OUTPUT
% MidImage: rotated midsagittal plane image
% SurfImage: rotated surface image
% xShift:
% yShift:
% 
%
% Author (contact): Clarissa Scholes (DePace Lab)
% Created: 2017-09-11
% Last Updated: 2020-07-27
% Documented by: Meghan Turner (meghan_turner@berkeley.edu)
%
% Functionalized from code originally in FindAPAxisFullEmbryo, by Meghan
% Turner (meghan_turner@berkeley) on 2020-07-27


%% CS20170911 This is really long-winded atm! Need to simplify.

%Find and open the mid-saggistal .nd files in the FullEmbryo folder
if find(~cellfun('isempty',strfind({dirFullEmbryo.name},'ANT_mid')))>0
    SurfFileIndexAnt = find(~cellfun('isempty',strfind({dirFullEmbryo.name},'ANT_mid')));
    SurfFileIndexPost = find(~cellfun('isempty',strfind({dirFullEmbryo.name},'POST_mid')));

    SurfImageAnt = bfopen([fullEmbryoPath,dirFullEmbryo(SurfFileIndexAnt).name]);
    SurfAntMeta = SurfImageAnt{4};

    SurfImagePost = bfopen([fullEmbryoPath,dirFullEmbryo(SurfFileIndexPost).name]);
    SurfPostMeta = SurfImagePost{4};

else
    error('FullEmbryo images mislabeled? Should be of form ANT_mid, POST_mid, ANT_surf, POST_surf')
end

%Get some useful parameters
try
    NSlices = SurfAntMeta.getPixelsSizeZ(0).getValue();
    NChannels = SurfAntMeta.getPixelsSizeC(0).getValue();
    PizelSize = SurfAntMeta.getPixelsPhysicalSizeX(0).getValue();
catch
    NSlices = str2double(SurfAntMeta.getPixelsSizeZ(0));
    NChannels = str2double(SurfAntMeta.getPixelsSizeC(0));
    PizelSize = str2double(SurfAntMeta.getPixelsPhysicalSizeX(0).value);%Size in microns of one pixel. NB for some of the OME metadata you need to add .value on end to get it out
end

%Figure out which channel to use (His channel)
HisChannel=find(~cellfun(@isempty,strfind(lower({Channel1{1},Channel2{1}}),'mcherry'))|...
    ~cellfun(@isempty,strfind(lower({Channel1{1},Channel2{1}}),'his')));

%Get the stack of images in the histone channel for ANT and POST images and do a maximum projection
temp = SurfImageAnt{1}(:,2);
IndexAntPost = find(~cellfun('isempty',strfind({temp{1:end}},['C=', num2str(HisChannel), '/', num2str(NChannels)])));
AntImStack = SurfImageAnt{1}(IndexAntPost, 1);
PostImStack = SurfImagePost{1}(IndexAntPost, 1);

%Max Project His channel for ANT and POST images
for i = 1:NSlices
    AntImStack2(:,:,i) = AntImStack{i};
    PostImStack2(:,:,i) = PostImStack{i};
end
AntMaxProj = max(AntImStack2, [], 3); 
PostMaxProj = max(PostImStack2, [], 3);

%Write these to files because EmbryoStitchNoMargin (which will stitch them together) wants it that way
imwrite(AntMaxProj, [fullEmbryoPath,'AntMaxProj.tif']); 
imwrite(PostMaxProj, [fullEmbryoPath,'PostMaxProj.tif']); 

%While at it, EmbryoStitchNoMargin also needs a flatfield. I don't take
%one at this magnification, so make it up
FFImage = ones(512,672); %make this automatically find the size later
imwrite(FFImage, [fullEmbryoPath,'FF20x.tif']);
Margin = 0;

%Define inputs to EmbryoStitchNoMargin, which will stitch the ANT and POST files together and record the xShift and yShift needed to do so
FFFile=[fullEmbryoPath,'FF20x.tif'];  
LeftImageFile = [fullEmbryoPath,'AntMaxProj.tif'];
RightImageFile = [fullEmbryoPath,'PostMaxProj.tif']; 

if FlipAP==0    %If FlipAP not specified in input, assume Anterior image on Left and Posterior on Right. 
    [APImage,xShift,yShift]=EmbryoStitchNoMargin(LeftImageFile,RightImageFile, FFFile, [], [], Margin);
else          %Otherwise assume the opposite. 
    [APImage,xShift,yShift]=EmbryoStitchNoMargin(RightImageFile, LeftImageFile, FFFile, [], [], Margin);
end

%Rename for use later in code
MidImage = APImage; 
%Display image to check it's stitched correctly
imshow(imadjust(mat2gray(MidImage,[0 65535])),'DisplayRange',[],'InitialMagnification',100);

%While we're at it, AddParticlePosition also needs the surface embryo,
%so generate a stiched together max projection of the ANT and POST
%surface images as well. 
if find(~cellfun('isempty',strfind({dirFullEmbryo.name},'ANT_surf')))>0
    SurfFileIndexAnt = find(~cellfun('isempty',strfind({dirFullEmbryo.name},'ANT_surf.nd')));
    SurfFileIndexPost = find(~cellfun('isempty',strfind({dirFullEmbryo.name},'POST_surf.nd')));
    SurfImageAnt = bfopen([fullEmbryoPath,dirFullEmbryo(SurfFileIndexAnt).name]);
    SurfAntMeta = SurfImageAnt{4};
    SurfImagePost = bfopen([fullEmbryoPath,dirFullEmbryo(SurfFileIndexPost).name]);
    SurfPostMeta = SurfImagePost{4};   
else
    error('FullEmbryo images mislabeled? Should be of form ANT_mid, POST_mid, ANT_surf, POST_surf')
end

%Get some useful parameters
try
    NSlices = SurfAntMeta.getPixelsSizeZ(0).getValue();
    NChannels = SurfAntMeta.getPixelsSizeC(0).getValue();
    PizelSize = SurfAntMeta.getPixelsPhysicalSizeX(0).getValue();
catch
    NSlices = str2double(SurfAntMeta.getPixelsSizeZ(0));
    NChannels = str2double(SurfAntMeta.getPixelsSizeC(0));
    PizelSize = str2double(SurfAntMeta.getPixelsPhysicalSizeX(0).value);%Size in microns of one pixel. NB for some of the OME metadata you need to add .value on end to get it out
end

%Get the stack of images in the histone channel for ANT and POST images and do a maximum projection
temp = SurfImageAnt{1}(:,2);
IndexAntPost = find(~cellfun('isempty',strfind({temp{1:end}},['C=', num2str(HisChannel), '/', num2str(NChannels)])));
AntImStack = SurfImageAnt{1}(IndexAntPost, 1);
PostImStack = SurfImagePost{1}(IndexAntPost, 1);

%Max Project His channel for ANT and POST images
for i = 1:NSlices
    AntImStack2(:,:,i) = AntImStack{i};
    PostImStack2(:,:,i) = PostImStack{i};
end
AntMaxProj = max(AntImStack2, [], 3); 
PostMaxProj = max(PostImStack2, [], 3);

%Write these to files because EmbryoStitchNoMargin (which will stitch them together) wants it that way
%The FlatField has already been written, above
imwrite(AntMaxProj, [fullEmbryoPath,'AntMaxProj.tif']); 
imwrite(PostMaxProj, [fullEmbryoPath,'PostMaxProj.tif']); 

%Define inputs to EmbryoStitchNoMargin, which will stitch the ANT and POST files together and record the xShift and yShift needed to do so
LeftImageFile = [fullEmbryoPath,'AntMaxProj.tif'];
RightImageFile = [fullEmbryoPath,'PostMaxProj.tif']; 

if FlipAP==0    %If FlipAP not specified in input, assume Anterior image on Left and Posterior on Right. 
    [APImage,xShift,yShift]=EmbryoStitchNoMargin(LeftImageFile,RightImageFile, FFFile, [], [], Margin);
else          %Otherwise assume the opposite. 
    [APImage,xShift,yShift]=EmbryoStitchNoMargin(RightImageFile, LeftImageFile, FFFile, [], [], Margin);
end

%Rename for use later in code
SurfImage = APImage; 
    