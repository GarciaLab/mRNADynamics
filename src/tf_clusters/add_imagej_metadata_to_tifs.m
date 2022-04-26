function add_imagej_metadata_to_tifs(preproc_data_folder, experiment_prefix)

% DESCRIPTION
% This function generates new versions of the PreProcessedData tiff stacks
% which will have FIJI-specific metadata associated with them. This change 
% makes the tiffs more widely compatible with things such as 
% Python's tifffile package, which expects to find FIJI metadata.
%
% ARGUMENTS
% preproc_data_folder: Path to your PreProcessedData folder
% experiment_prefix: Prefix for the specific dataset
%
% OPTIONS
% none
%
% OUTPUT
% Replaces all PreProcessedData tiffs with identical tiffs that now have
% FIJI metadata associated with them
%
% Author (contact): Meghan Turner (meghanaturner1@gmail.com)
% Created: 04/19/2022
% Last Updated:
%% 

curr_dir = [preproc_data_folder, filesep, experiment_prefix];
all_tifs = dir([curr_dir, filesep, '*_ch0*tif']);

wbar = waitbar(0, ['Adding ImageJ metadata to ' experiment_prefix ' PreProc tiffs']);
for i = 1:length(all_tifs)
    waitbar(i/length(all_tifs), wbar);

    curr_tif = tiffreadVolume([curr_dir, filesep, all_tifs(i).name]);
    im_stack = curr_tif; 
    
    % Make a new Tiff object, name it the same as the old tif file
    t = Tiff([curr_dir, filesep, all_tifs(i).name],'w');
    
    fiji_descr = ['ImageJ=1.52u', newline, ...
                  'images=', num2str(size(im_stack,3)), newline,... 
                  'slices=', num2str(size(im_stack,3)), newline,...
                  'loop=false', newline,...  
                  'min=0.0', newline,...      
                  'max=', num2str(max(im_stack(:)))];  % 16-bit image
    %             'channels=' num2str(1) newline...
    %             'frames=' num2str(1) newline... 
    %             'hyperstack=true' newline...
    %             'mode=grayscale' newline... 

    % Add the ImageJ/FIJI metadats to the Tiff tag structure
    tagstruct.ImageDescription = fiji_descr;
    % 7 other mandatory tags per the Tiff class documentation:
    tagstruct.ImageLength = size(im_stack,1);
    tagstruct.ImageWidth = size(im_stack,2);
    tagstruct.Photometric = Tiff.Photometric.MinIsBlack;
    tagstruct.BitsPerSample = 16;
    tagstruct.SamplesPerPixel = 1;
    tagstruct.Compression = Tiff.Compression.LZW; % DO NOT use None - it circularly permutes the pixels in xy space
    tagstruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
    % optional tag:
    tagstruct.SampleFormat = Tiff.SampleFormat.UInt;
    
    % Now write the image stack to the new Tiff object
    for slice = 1:size(im_stack,3)
        t.setTag(tagstruct)
        t.write(im2uint16(im_stack(:,:,slice)));
        t.writeDirectory(); % saves a new page in the tiff file
    end
    t.close()
end
close(wbar)
disp('Done!')