function segment_track_spots(Prefix, thresh, show_status, save_status, ...
    track_status, num_frames)

% show_status takes 0 or 1 depending on whether you want to display plots
% and images. When using show_status = 1, change "parfor" by "for" in the
% loop that goes over all spots

% save_status takes 0 or 1 depending on whether you want to save your files

% track_status takes 0 or 1 depending on whether you want to make a time
% tracking in addition to the segmentation.

% num_frames for debugging should be kept at 5-20

% thresh should be kept at ~90-200 for lattice data, and at ~5-10 for
% confocal data.

%% 

tic;

[SourcePath,FISHPath,DropboxFolder,MS2CodePath,PreProcPath]=...
    DetermineLocalFolders(Prefix);

%Here I want to load each of the tifs and iterate through them, generating
%dog images in each iteration

% Load data generated with ExportDataForFISH


load([DropboxFolder,filesep,Prefix,filesep,'FrameInfo.mat']);
im_stack = {};

if num_frames == 0
    num_frames = length(FrameInfo);
end

for current_frame = 1:num_frames
    for current_z = 1:FrameInfo(1).NumberSlices
        im_stack{current_frame,current_z} = imread([PreProcPath,filesep,Prefix,filesep,Prefix,'_',iIndex(current_frame,3),'_z',iIndex(current_z,2),'.tif']);
    end
end

OutputFolder1 = [FISHPath,filesep,Prefix,'_',filesep,'MYCODEdogsMYCODE'];
OutputFolder2 = [FISHPath,filesep,Prefix,'_',filesep,'MYCODEsegsMYCODE'];

mkdir(OutputFolder1)
mkdir(OutputFolder2)
%%
%DoG Stuff
%filterSize >> sigma 2 > sigma 1. these values should be good for a first
%pass.

pixelSize = FrameInfo(1).PixelSize*1000; %nm

sigma1 = 150 / pixelSize;
sigma2 = 250 / pixelSize;
filterSize = round(1500 / pixelSize); 
neighb = round(500 / pixelSize); %This should work for a first pass and shouldn't fail on sisters.
thr = thresh; 
dog_stack  = {}; 
all_frames = {}; 
for current_frame = 1:num_frames-1 
    for current_z = 1:size(im_stack,2) %z-slices
        im = im_stack{current_frame,current_z};
        %filterSize >> sigma 2 > sigma 1. these values should be good for a first pass.
        dog_stack{current_frame,current_z} = conv2(single(im), single(fspecial('gaussian',filterSize, sigma1) - fspecial('gaussian',filterSize, sigma2)),'same');
        if save_status
            dog_name = ['DOG_',Prefix,'_',iIndex(current_frame,3),'_z',iIndex(current_z,2),'.tif'];
            imwrite(uint16(dog_stack{current_frame,current_z}), [OutputFolder1,filesep,dog_name])
        end
        dog = dog_stack{current_frame,current_z}(10:end-10, 10:end-10);
        im = im(10:end-10, 10:end-10);
        if show_status
            figure(1)
            imshow(im,[]);
        end
        thrim = dog > thr;
        [im_label, n_spots] = bwlabel(thrim); 
        temp_particles = {};
%         rad = 500/pixelSize; %500nm is roughly the size of a sister chromatid diffraction limited spot.
        rad = 800/pixelSize;
        temp_frames = {};
        if n_spots ~= 0
            parfor k = 1:n_spots
                temp_particles(k) = fit_single_spot(k, im, im_label, dog, ...
                    neighb, rad, pixelSize, show_status);
                if k == n_spots && save_status
                    seg_name = ['SEG_',Prefix,'_',iIndex(current_frame,3),'_z',iIndex(current_z,2),'.tif'];
                    saveas(gcf,[OutputFolder2,filesep,seg_name]);
                end
            end
            for k = 1:n_spots
                if ~isempty(temp_particles{k})
                    temp_frames = [temp_frames, temp_particles(k)];
                end
            end
            all_frames{current_frame,current_z} = temp_frames;
        end
    end
end

%%

clear im_stack
clear dog_stack
n = 1;
nframes = size(all_frames,1);
nz = size(all_frames,2);
for i = 1:nframes 
    for j = 1:nz 
         for spot = 1:length(all_frames{i,j}) %spots within particular image
             if ~isempty(all_frames{i,j}{spot})
                 Particles(n).CentralIntensity(1) = cell2mat(all_frames{i,j}{spot}(12));
                 Particles(n).FixedAreaIntensity(1) = cell2mat(all_frames{i,j}{spot}(1));
                 Particles(n).GaussianIntensity(1) = cell2mat(all_frames{i,j}{spot}(11));
                 Particles(n).DOGIntensity(1) = cell2mat(all_frames{i,j}{spot}(13));
                 Particles(n).xFit(1) = cell2mat(all_frames{i,j}{spot}(2));
                 Particles(n).yFit(1) = cell2mat(all_frames{i,j}{spot}(3));
                 Particles(n).xDoG(1) = cell2mat(all_frames{i,j}{spot}(10));
                 Particles(n).yDoG(1) = cell2mat(all_frames{i,j}{spot}(9));
                 Particles(n).Offset(1) = cell2mat(all_frames{i,j}{spot}(4));
%                  Particles(n).Sister(1) = all_frames{i,j}{spot}(5);
                 Particles(n).Snippet{1} = cell2mat(all_frames{i,j}{spot}(5));
                 Particles(n).Area{1} = cell2mat(all_frames{i,j}{spot}(6));
                 Particles(n).xFitWidth{1} = cell2mat(all_frames{i,j}{spot}(7));
                 Particles(n).yFitWidth{1} = cell2mat(all_frames{i,j}{spot}(8));
                 Particles(n).z(1) = j;
                 Particles(n).frame(1) = i;
                 Particles(n).r = 0;
                 n = n + 1;
             end
         end
    end
end

fields = fieldnames(Particles);

% z tracking
changes = 1;
while changes ~= 0
    changes = 0;
    i = 1;
    for n = 1:nframes   
        i = i + length(Particles([Particles.frame] == (n - 1) ));
        for j = i:i+length(Particles([Particles.frame] == n)) - 1
            for k = j+1:i+length(Particles([Particles.frame] == n)) - 1
                dist = sqrt( (Particles(j).xFit(end) - Particles(k).xFit(end))^2 + (Particles(j).yFit(end) - Particles(k).yFit(end))^2); 
                if dist < neighb && Particles(j).z(end) ~= Particles(k).z(end)
                    for m = 1:numel(fields)-2 %do not include fields 'r' or 'frame'
                        Particles(j).(fields{m}) = [Particles(j).(fields{m}), Particles(k).(fields{m})];
                    end
                    Particles(k).r = 1;
                    changes = changes + 1;
                end
            end
        end
    end
    Particles = Particles([Particles.r]~=1);
end


%pick the brightest z-slice
for i = 1:length(Particles)
    [~, max_index] = max(Particles(i).FixedAreaIntensity);
    if track_status
        for j = 1:numel(fields)-2 %do not include fields 'r' or 'frame'
            Particles(i).(fields{j}) = Particles(i).(fields{j})(max_index);
        end
    else 
        Particles(i).brightestZ = Particles(i).z(max_index);
    end
end

%Create a final Spots structure to be fed into TrackmRNADynamics

Spots = [];
for i = 1:Particles(end).frame
    frames = find([Particles.frame]==i);
    if ~isempty(frames)
        for j = frames(1):frames(end)
            Spots(i).Fits(j-frames(1)+1) = Particles(j);
        end
    end
end

%time tracking

if track_status
    Particles = track_spots(Particles, neighb);
end

%Save and plot

mkdir([DropboxFolder,filesep,Prefix,filesep,'mycode']);
save([DropboxFolder,filesep,Prefix,filesep,'mycode',filesep,'Spots.mat'], 'Spots');
for i = 1:length(Particles)
    if length(Particles(i).frame) > 50
        plot(Particles(i).frame, Particles(i).FixedAreaIntensity)
        i
        hold on
    end
end

t = toc