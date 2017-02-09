function segmentSpotsML(Prefix,Threshold,varargin)

%Parameters:
%Prefix: Prefix of the data set to analyze
%Threshold: Threshold to be used. Should be kept at ~90-200 for lattice
%           light-sheet data, and at ~5-10 for confocal data (Leica SP8).
%           If left empty, then the code just generates the DoG files.

%Options:
%'displayFigures':   If you want to display plots and images.
%                
%'TrackSpots':   Do you want to use this code to track the particles instead
%                of using TrackmRNADynamics? 
%'Frames',N:     Run the code from frame 1 to frame N. Defaults to all
%                frames. It's suggested to run 5-20 frames for debugging.
%'NoShadows':    Spots without valid (peaked) z-profiles are normally
%                discarded. This option overrides that.
%               

%Default options
displayFigures=0;
TrackSpots=0;
num_frames=0;
Shadows = 1;

for i=1:length(varargin)
    if strcmp(varargin{i},'displayFigures')
        displayFigures=1;
    elseif strcmp(varargin{i},'TrackSpots')
        TrackSpots=1;
    elseif strcmpi(varargin{i}, 'NoShadows')
        Shadows = 0;
    elseif strcmp(varargin{i},'Frames')
        if ~isnumeric(varargin{i+1})
            error('Wrong input parameters. After ''Frames'' you should input the number of frames')
        else
            num_frames=varargin{i+1};
        end
    end
end

%If no threshold was specified, then just generate the DoG images
try
    if isempty(Threshold)
        just_dog=1;
    else
        just_dog=0;
    end
catch
    error('Please pass the argument "[]" to generate DoG images')
end
%%
tic;

[~,~,~,~,~,~,~,ExperimentType, Channel1, Channel2,~] =...
    readMovieDatabase(Prefix);

[SourcePath,FISHPath,DropboxFolder,MS2CodePath,PreProcPath]=...
    DetermineLocalFolders(Prefix);

load([DropboxFolder,filesep,Prefix,filesep,'FrameInfo.mat']);
microscope = FrameInfo(1).FileMode;
zSize = FrameInfo(1).NumberSlices + 2;
if num_frames == 0
    num_frames = length(FrameInfo);
end
OutputFolder1=[FISHPath,filesep,Prefix,'_',filesep,'dogs'];
mkdir(OutputFolder1)

%Load flat-field
doFF = 1;
try
    ffim = imread([PreProcPath, filesep, Prefix, filesep,Prefix,'_FF.tif']);
    ffim = CPsmooth(ffim,'Gaussian Filter',256,0);
    ffim = double(ffim/max(max(ffim)));
catch
    warning('Will not apply flat field correction');
    doFF = 0;
end

clear rawdir;
%%

%The spot finding algorithm first segments the image into regions that are
%above the threshold. Then, it finds global maxima within these regions by searching in a region "neighborhood"
%within the regions. 


pixelSize = FrameInfo(1).PixelSize*1000; %nm
neighborhood = round(1300 / pixelSize); %nm
snippet_size = 2*(floor(1300/(2*pixelSize))) + 1; % nm. note that this is forced to be odd

           
all_frames = cell(num_frames, zSize);
close all force;
if just_dog
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Generate difference of Gaussian images if no threshold was given.
    %Initialize Difference of Gaussian filter parameters. filterSize >> sigma2
    %> sigma1
    sigma1 = pixelSize / pixelSize; %width of narrower Gaussian
    sigma2 = 42000 / pixelSize; % width of wider Gaussian
    filterSize = round(2000 / pixelSize); %size of square to be convolved with microscopy images
    zim = [];
    try
    %this is just some function that can only be called if IJM is set up
    IJM.getIdentifier() 
    catch
        addpath('e:/Fiji.app/scripts') % Update for your ImageJ installation
        ImageJ                         % Initialize IJM and MIJ
    end
    evalin('base', 'IJM')
    evalin('base', 'MIJ')
for current_frame = 1:num_frames
    for i = 1:zSize
        zim(:,:,i) = double(imread([PreProcPath,filesep,Prefix, filesep, Prefix,'_',iIndex(current_frame,3),'_z',iIndex(i,2),'.tif']));
    end
    mkdir([PreProcPath,filesep,Prefix,filesep,'stacks']);
    name = [PreProcPath,filesep,Prefix,filesep,'stacks', filesep, iIndex(current_frame,3),'.tif'];
    imwrite(uint16(zim(:,:,1)), name);
    for k = 2:size(zim,3)
        imwrite(uint16(zim(:,:,k)), name, 'WriteMode', 'append');
    end
    mij.run('Trainable Weka Segmentation 3D', ['open=',name]);
    pause(20);
    trainableSegmentation.Weka_Segmentation.loadClassifier('C:\\Users\\ArmandoReimer\\Desktop\\weca\\vasa_first_movie_3dmachine\\classifier.model');
    trainableSegmentation.Weka_Segmentation.getProbability();
    ijm.getDatasetAs('probmaps')
    p = evalin('base', 'probmaps');
    p2 = [];
    for m = 1:2:46
        p2(:,:,ceil(m/2)) =  p(:,:,m); %the even images in the original array are negatives of the odds
    end
    p2 = permute(p2, [2 1 3]);
    p2 = p2*10000;
    for i = 1:size(p2, 3)      
        p_name = ['prob',Prefix,'_',iIndex(current_frame,3),'_z',iIndex(i,2),'.tif'];
        imwrite(uint16(p2(:,:,i)), [OutputFolder1,filesep,p_name])
    end
    MIJ.run('Close All');
end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Segment transcriptional loci
else
    h=waitbar(0,'Segmenting spots');
    for current_frame = 1:num_frames
        w = waitbar(current_frame/num_frames,h);
        set(w,'units', 'normalized', 'position',[0.4, .15, .25,.1]);
        for i = 1:zSize   
            if strcmpi(ExperimentType, 'inputoutput')
                im = double(imread([PreProcPath,filesep,Prefix,filesep,Prefix,'_',iIndex(current_frame,3),'_z',iIndex(i,2),'_ch02','.tif']));
            else
                im = double(imread([PreProcPath,filesep,Prefix,filesep,Prefix,'_',iIndex(current_frame,3),'_z',iIndex(i,2),'.tif']));
            end
            dog = double(imread([OutputFolder1, filesep,'prob',Prefix,'_',iIndex(current_frame,3),'_z',iIndex(i,2),'.tif']));
%             if displayFigures
%                 fig = figure(1);
%                 imshow(dog,[]);
%             else
%                 fig=[];
%             end
            %apply flatfield correction
            if doFF && sum(size(im)==size(ffim))
                im = im./ffim;
            else
                warning('Not applying flat-field correction to MS2 channel')
            end
            %
            im_thresh = dog >= Threshold;
%             se = strel('square', 4);
%             im_thresh = imdilate(im_thresh, se); %thresholding from this classified probability map can produce non-contiguous, spurious spots. This fixes that and hopefully does not combine real spots from different nuclei
            im_thresh = imgaussfilt(uint8(im_thresh)*255, 1);
            im_thresh = im_thresh>0;
            [im_label, n_spots] = bwlabel(im_thresh); 
              
            if displayFigures
                fig = figure(1);
                imshow(im_thresh,[]);
            else 
                fig = [];
            end
                           
            temp_frames = {};
            temp_particles = cell(1, n_spots);

            if n_spots ~= 0
                if ~displayFigures
                    parfor k = 1:n_spots
                        try
                            temp_particles(k) = identifySingleSpot(k, im, im_label, dog, ...
                                neighborhood, snippet_size, pixelSize, displayFigures, fig, microscope, 0);
                        catch
                        end
                    end
                else
                    for k = 1:n_spots
                            temp_particles(k) = identifySingleSpot(k, im, im_label, dog, ...
                                neighborhood, snippet_size, pixelSize, displayFigures, fig, microscope, 0);
                    end
                end
                for k = 1:n_spots
                    if ~isempty(temp_particles{k})
                        temp_frames = [temp_frames, temp_particles(k)];
                    end
                end
                all_frames{current_frame,i} = temp_frames;
            end
        end
    end
    close(h)
    close all force;
end

%%
%Create a useful structure that can be fed into pipeline
if ~just_dog 
    n = 1;
    h=waitbar(0,'Saving particle information');
    for i = 1:num_frames  
        waitbar(i/num_frames,h)
        for j = 1:zSize 
             for spot = 1:length(all_frames{i,j}) %spots within particular image
                 if ~isempty(all_frames{i,j}{spot})
                     Particles(n).FixedAreaIntensity(1) = cell2mat(all_frames{i,j}{spot}(1));
                     Particles(n).xFit(1) = cell2mat(all_frames{i,j}{spot}(2));
                     Particles(n).yFit(1) = cell2mat(all_frames{i,j}{spot}(3));
                     Particles(n).Offset(1) = cell2mat(all_frames{i,j}{spot}(4));
                     Particles(n).Snippet{1} = cell2mat(all_frames{i,j}{spot}(5));
                     Particles(n).Area{1} = cell2mat(all_frames{i,j}{spot}(6));
                     Particles(n).xFitWidth{1} = cell2mat(all_frames{i,j}{spot}(7));
                     Particles(n).yFitWidth{1} = cell2mat(all_frames{i,j}{spot}(8));
                     Particles(n).yDoG(1) = cell2mat(all_frames{i,j}{spot}(9));
                     Particles(n).xDoG(1) = cell2mat(all_frames{i,j}{spot}(10));
                     Particles(n).GaussianIntensity(1) = cell2mat(all_frames{i,j}{spot}(11));                     
                     Particles(n).CentralIntensity(1) = cell2mat(all_frames{i,j}{spot}(12));
                     Particles(n).DOGIntensity(1) = cell2mat(all_frames{i,j}{spot}(13));
                     Particles(n).snippet_mask{1} = cell2mat(all_frames{i,j}{spot}(14));
                     Particles(n).SisterDistance(1) = cell2mat(all_frames{i,j}{spot}(17));
                     Particles(n).ConfidenceIntervals{1} = cell2mat(all_frames{i,j}{spot}(19));          
                     Particles(n).gaussSpot{1} = cell2mat(all_frames{i,j}{spot}(20));
                     raw = all_frames{i,j}{spot}(21);
                     Particles(n).rawSpot{1} = raw{1};
                     Particles(n).z(1) = j;
                     Particles(n).discardThis = 0;
                     Particles(n).frame(1) = i;
                     Particles(n).r = 0;
                     n = n + 1;
                 end
             end
        end
    end
    close(h)
    fields = fieldnames(Particles);
    
    %z-tracking
    changes = 1;
    while changes ~= 0
        changes = 0;
        i = 1; 
        h=waitbar(0,'Finding z-columns');
        neighborhood = 1700 / pixelSize;
        for n = 1:num_frames  
            waitbar(n/num_frames,h)
            l = length(Particles([Particles.frame] == n));
            i = i + length(Particles([Particles.frame] == (n - 1) ));
            for j = i:i+l-1
                for k = j+1:i+l-1
                    dist = sqrt( (Particles(j).xFit(end) - Particles(k).xFit(end))^2 + (Particles(j).yFit(end) - Particles(k).yFit(end))^2); 
                    if dist < neighborhood && Particles(j).z(end) ~= Particles(k).z(end)
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
        close(h)
    end
    falsePositives = 0;
    %pick the brightest z-slice
    for i = 1:length(Particles)
        [~, max_index] = max(Particles(i).CentralIntensity);
        if TrackSpots
            for j = 1:numel(fields)-2 %do not include fields 'r' or 'frame'
                Particles(i).(fields{j}) = Particles(i).(fields{j})(max_index);
            end
        else
            Particles(i).brightestZ = Particles(i).z(max_index);
            if Shadows && Particles(i).brightestZ == min(Particles(i).z)...
               || Particles(i).brightestZ == max(Particles(i).z)
                Particles(i).discardThis = 1;
                Particles(i).noIntensityAnalysis = 1;
                falsePositives = falsePositives + 1;
            end  
        end
    end
 
    %Create a final Spots structure to be fed into TrackmRNADynamics
    Spots = [];            
    fields = fieldnames(Particles);
    num_fields = length(fields);
    for i = 1:num_frames
        frames = find([Particles.frame]==i);
        if ~isempty(frames)
            for j = frames(1):frames(end)
                if ~Particles(j).discardThis
                    Spots(i).Fits(j-frames(1)+1) = Particles(j);
                end
                %Sometimes, all spots are discarded in a frame. In that
                %case, create an empty Spots entry in that frame.
                if length(Spots)<i
                    for l = 1:num_fields
                        Spots(i).Fits.(fields{l}) = [];
                    end
                end
            end
        else 
            for l = 1:num_fields
                Spots(i).Fits.(fields{l}) = [];
            end
        end
    end
    
    %Clean up Spots to remove empty rows
    Spots2 = struct('Fits', []);
    for i = 1:length(Spots)
        Spots2(i).Fits = [];
        for j = 1:length(Spots(i).Fits)
            if j~=1
                if ~isempty(Spots(i).Fits(j).z)...
                        && ~isequal(Spots(i).Fits(j).CentralIntensity,Spots(i).Fits(j-1).CentralIntensity)
                        Spots2(i).Fits = [Spots2(i).Fits, Spots(i).Fits(j)];                    
                end
            else
                if ~isempty(Spots(i).Fits(j).z)                        
                    Spots2(i).Fits = [Spots2(i).Fits, Spots(i).Fits(j)];
                end
            end
        end
    end
    Spots = Spots2;
 
    %AR 7/10/16: Optional time tracking using track_spots script. Also
    %makes some potentially useful plots. This was originally here to have
    %a single, fully integrated script before this segmentation was worked
    %into the rest of the pipeline.
    neighborhood = 3000 / pixelSize;
    if TrackSpots
        Particles = track_spots(Particles, neighborhood);
        save([DropboxFolder,filesep,Prefix,filesep,'Particles_SS.mat'], 'Particles');
    end

    mkdir([DropboxFolder,filesep,Prefix]);
    save([DropboxFolder,filesep,Prefix,filesep,'Spots.mat'], 'Spots');    

end

t = toc;
display(['Elapsed time: ',num2str(t/60),' min'])

if ~just_dog
    detectedCircles = 0;
    detectedBalls = 0;
    for i = 1:length(Spots)
        for j = 1:length(Spots(i).Fits)
            detectedCircles = detectedCircles + length(Spots(i).Fits(j).z);
            detectedBalls = detectedBalls + 1;
        end
    end
    display(['Detected spots: ',num2str(detectedCircles)])
    if exist([DropboxFolder,filesep,Prefix,filesep,'log.mat'])
        load([DropboxFolder,filesep,Prefix,filesep,'log.mat'])
        log(end+1).Date = date;
        log(end).runTime = t/60; %min
        log(end).falsePositives = falsePositives;
        log(end).totalCircles = detectedCircles;
        log(end).totalBalls = detectedBalls;
        log(end).avgZSize = detectedCircles/detectedBalls;
        log(end).Threshold = Threshold;
    else
        log = struct();
        log(1).Date = date;
        log(1).runTime = t/60; %min
        log(1).falsePositives = falsePositives;
        log(1).totalCircles = detectedCircles;
        log(1).totalBalls = detectedBalls;
        log(1).avgZSize = detectedCircles/detectedBalls;
        log(1).Threshold = Threshold;
    end
    save([DropboxFolder,filesep,Prefix,filesep,'log.mat'], 'log');
end
try
    poolobj = gcp('nocreate');
    delete(poolobj);
catch
end