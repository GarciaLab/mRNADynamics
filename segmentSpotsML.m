function segmentSpotsML(Prefix,Threshold,varargin)
% segmentSpotsML(Prefix, Threshold, [Options])
%
% DESCRIPTION
% Identify and segment individual transcription Spots. 
%
% ARGUMENTS
% Prefix: Prefix of the data set to analyze
% Threshold: Threshold to be used. Should be kept at 5000 always.
%           If left empty, then the code just generates the DoG files.
% [Options]: See below.
%
% OPTIONS
% 'displayFigures':   If you want to display plots and images.
% 'Tifs':         When running this script without a threshold to generate
%                probability maps, use this option to instead only generate
%                the TIF stacks necessary for doing Weka classification. 
%                Recommended to run this before making a new classifier.
%                
% 'TrackSpots':   Do you want to use this code to track the particles instead
%                of using TrackmRNADynamics?
% 'InitialFrame', N: Run the code from frame N to last frame. Defaults to first
%                frame.
% 'LastFrame', M:     Run the code from initial frame to frame M. Defaults to all
%                frames. It's suggested to run 5-20 frames for debugging.
% 'Shadows':    	 This option should be followed by 0, 1 or 2. This
%                specifies the number of requisite z-planes above and/or below the
%                brightest plane for a spot to have to pass quality control. 
%               
% OUTPUT
% 'Spots':  A structure array with a list of detected transcriptional loci
% in each frame and their properties.
%
% Author (contact): Armando Reimer (areimer@berkeley.edu)
% Created: 01/01/2016
% Last Updated: 12/31/2017
%
% Documented by: Armando Reimer (areimer@berkeley.edu)


%Default options
displayFigures=0;
TrackSpots=0;
num_frames=0;
num_shadows = 2;
initial_frame = 1;
just_tifs = 0;
    
for i=1:length(varargin)
    if strcmp(varargin{i},'displayFigures')
        displayFigures=1;
    elseif strcmp(varargin{i}, 'Tifs')
        just_tifs = 1;
    elseif strcmp(varargin{i},'TrackSpots')
        TrackSpots=1;
    elseif strcmp(varargin{i},'LastFrame')
        if ~isnumeric(varargin{i+1})
            error('Wrong input parameters. After ''Frames'' you should input the number of frames')
        else
            num_frames=varargin{i+1};
        end
    elseif strcmp(varargin{i},'Shadows')
        if ~isnumeric(varargin{i+1}) || varargin{i+1} > 2
            error('Wrong input parameters. After ''Shadows'' you should input number of shadows (0, 1 or 2)')
        else
            num_shadows=varargin{i+1};
        end
    elseif strcmp(varargin{i}, 'InitialFrame')
         if ~isnumeric(varargin{i+1}) || varargin{i+1} < 1
            error('Wrong input parameter for initial frame.')
        else
            initial_frame=varargin{i+1};
        end
    else
        if ~isnumeric(varargin{i})
            error('Input parameters not recognized. Check spelling and case.')
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

maxWorkers = 6;
try
    parpool(maxWorkers);  % 6 is the number of cores the Garcia lab server can reasonably handle per user at present.
catch
    %parpool throws an error if there's a pool already running. 
end
    
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

nCh = 1;
if strcmpi(ExperimentType, '2spot2color')
    nCh = 2;
end

if strcmpi(ExperimentType,'inputoutput')
    if  contains(Channel2,'mcp', 'IgnoreCase', true) ||...
            contains(Channel2,'pcp','IgnoreCase',true)
        coatChannel=2;
    elseif  contains(Channel1,'mcp', 'IgnoreCase', true) ||...
            contains(Channel1,'pcp','IgnoreCase',true)
        coatChannel=1;
    else
        error('No MCP or PCP channel detected. Check MovieDatabase.XLSX')
    end
end
            
%Load and apply flat-field correction
doFF = 1;
try
    ffim = imread([PreProcPath, filesep, Prefix, filesep,Prefix,'_FF.tif']);
    ffim = CPsmooth(ffim,'Gaussian Filter',256,0); %large feature scale-space representation
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
%Generate probability maps of likely transcriptional loci

stacksPath = [PreProcPath,filesep,Prefix,filesep,'stacks', ];
mkdir(stacksPath);
rawStackArray = [];
if ~just_tifs    
    [classifierPathCh1,classifierFolder]=uigetfile([MS2CodePath, filesep, 'classifiers', filesep, '*.model']);
    if nCh==2
      [classifierPathCh2,~]=uigetfile([MS2CodePath, filesep, 'classifiers', filesep, '*.model']); 
    end
    evalin('base', 'clear probmaps');
    version -java;
    javaver = ans;
    if ~strcmp('Java 1.8.0_66-b18 with Oracle Corporation Java HotSpot(TM) 64-Bit Server VM mixed mode', javaver)
        error('Java version incorrect. Re-run InstallmRNADynamics or check environment variables')
    end
    heapSize = java.lang.Runtime.getRuntime.maxMemory;
    if heapSize<1E10 
        error('Please increase your Java heap memory allocation to at least 10GB (Home -> Preferences -> General -> Java Heap Memory.');
    end
    
    try
        %this is just some function that can only be called if IJM is set up
        IJM.getIdentifier() 
    catch
        addpath([MS2CodePath,filesep,'Fiji.app',filesep,'scripts'])
        ImageJ               % Initialize IJM and MIJ
    end
   
    ijm = evalin('base', 'IJM');
    mij = evalin('base', 'MIJ');
    zSize2 = zSize*2;
    h=waitbar(0,'Running Weka Classifier');
else 
    h = waitbar(0, 'Making .tif stacks for Weka classification');
end
%Make requisite TIF stacks for classification
for q = 1:nCh

    if strcmpi(ExperimentType,'inputoutput')               
        nameSuffix= ['_ch',iIndex(coatChannel,2)];
    else
        nameSuffix= ['_ch',iIndex(q,2)];
    end
    
    for current_frame = initial_frame:num_frames
        w = waitbar(current_frame/num_frames);
        set(w,'units', 'normalized', 'position',[0.4, .15, .25,.1]);
        for i = 1:zSize
            rawStackArray(:,:,i) = imread([PreProcPath,filesep,Prefix, filesep, Prefix,'_',iIndex(current_frame,3),'_z',iIndex(i,2),nameSuffix,'.tif']);    
        end
        rawStackName = [stacksPath, filesep, iIndex(current_frame,3), nameSuffix,'.tif'];
        %Don't write new stacks if they're already made.
        if length(dir([stacksPath, filesep, '*.tif'])) ~= num_frames
            imwrite(uint16(rawStackArray(:,:,1)), rawStackName);
            for k = 2:size(rawStackArray,3)
                imwrite(uint16(rawStackArray(:,:,k)), rawStackName, 'WriteMode', 'append');
            end
        end
        %Do the classification with Weka in Fiji
        if ~just_tifs
            mij.run('Trainable Weka Segmentation 3D', ['open=',rawStackName]);
            pause(10);
            if q==1 || strcmpi(ExperimentType,'inputoutput')
                trainableSegmentation.Weka_Segmentation.loadClassifier([classifierFolder, classifierPathCh1]);
            elseif q==2
                trainableSegmentation.Weka_Segmentation.loadClassifier([classifierFolder, classifierPathCh2]);
            else
                error(['This pipeline does not support',...
                    'more than two spot channels. If you''re actually',...
                    'trying to segment 3 or more channels, talk to Armando to',...
                    'get this implemented. Otherwise you''ve reached an error.',...
                    'Check your data. This is probably not a bug in the code.']);
            end
            trainableSegmentation.Weka_Segmentation.getProbability();
            ijm.getDatasetAs('probmaps')
            pMapTemp = evalin('base', 'probmaps');
            pMap = [];
            for m = 1:2:zSize2
                pMap(:,:,ceil(m/2)) =  pMapTemp(:,:,m); %the even images in the original array are negatives of the odds
            end
            pMap = permute(pMap, [2 1 3]) * 10000; %multiplying so this can be cast to uint16
            for i = 1:size(pMap, 3)
                p_name = ['prob',Prefix,'_',iIndex(current_frame,3),'_z',iIndex(i,2),nameSuffix,'.tif'];
                imwrite(uint16(pMap(:,:,i)), [OutputFolder1,filesep,p_name])               
            end
            mij.run('Close All');
        end
    end
    close(w);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Segment transcriptional loci
else
    for q=1:nCh

        if strcmpi(ExperimentType,'inputoutput')            
            nameSuffix= ['_ch',iIndex(coatChannel,2)];
        else
            nameSuffix = ['_ch',iIndex(q,2)];
        end
        
        h=waitbar(0,'Segmenting Spots');
        for current_frame = initial_frame:num_frames
            w = waitbar(current_frame/(num_frames-initial_frame),h);
            set(w,'units', 'normalized', 'position',[0.4, .15, .25,.1]);
            for i = 1:zSize   
                im = double(imread([PreProcPath,filesep,Prefix,filesep,Prefix,'_',iIndex(current_frame,3),'_z',iIndex(i,2),nameSuffix,'.tif']));
                pMap = double(imread([OutputFolder1, filesep,'prob',Prefix,'_',iIndex(current_frame,3),'_z',iIndex(i,2),nameSuffix,'.tif']));
                if displayFigures
                    fig = figure(1);
                    imshow(im,[]);
                else
                    fig=[];
                end
                %apply flatfield correction
                if doFF && sum(size(im)==size(ffim))
                    im = im.*ffim;
                end
                %
                im_thresh = pMap >= Threshold(q);
                se = strel('square', 3);
                im_thresh = imdilate(im_thresh, se); %thresholding from this classified probability map can produce non-contiguous, spurious Spots{q}. This fixes that and hopefully does not combine real Spots{q} from different nuclei
                im_thresh = im_thresh>0;
                [im_label, n_Spots] = bwlabel(im_thresh); 
                centroids = regionprops(im_thresh, 'centroid');
    %               
    %             if displayFigures
    %                 fig = figure(1);
    %                 imshow(im_thresh,[]);
    %             else 
    %                 fig = [];
    %             end

                temp_frames = {};
                temp_particles = cell(1, n_Spots);

                if n_Spots ~= 0
                    if ~displayFigures                    
                        parfor k = 1:n_Spots
                            try
                                centroid = round(centroids(k).Centroid);
                                temp_particles(k) = identifySingleSpot(k, im, im_label, pMap, ...
                                    neighborhood, snippet_size, pixelSize, displayFigures, fig, microscope, 0, centroid, 'ML');
                            catch 
                            end
                        end
                    else
                        for k = 1:n_Spots
                            centroid = round(centroids(k).Centroid);    
                            temp_particles(k) = identifySingleSpot(k, im, im_label, pMap, ...
                                neighborhood, snippet_size, pixelSize, displayFigures, fig, microscope, 0, centroid, 'ML');
                        end
                    end
                    for k = 1:n_Spots
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
        for i = initial_frame:num_frames  
            waitbar(i/(num_frames-initial_frame),h)
            for j = 1:zSize 
                 for spot = 1:length(all_frames{i,j}) %Spots{q} within particular image
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
        
        try
            fields = fieldnames(Particles);
        catch 
            error('No spots found, probably. Did you want spots? Try something else.')
        end
        
        %z-tracking
        changes = 1;
        while changes ~= 0
            changes = 0;
            i = 1; 
            h=waitbar(0,'Finding z-columns');
            neighborhood = 1300 / pixelSize;
            for n = initial_frame:num_frames
                waitbar(n/(num_frames-initial_frame),h)
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
                if num_shadows == 1
                    if length(Particles(i).z) <= 1
                            Particles(i).discardThis = 1;
                            Particles(i).noIntensityAnalysis = 1;
                            falsePositives = falsePositives + 1;
                    elseif  Particles(i).brightestZ == Particles(i).z(end)
                        if Particles(i).z(max_index -1) ~= Particles(i).brightestZ-1
                            Particles(i).discardThis = 1;
                            Particles(i).noIntensityAnalysis = 1;
                            falsePositives = falsePositives + 1;
                        end
                    elseif Particles(i).brightestZ == Particles(i).z(1)
                        if Particles(i).z(max_index+1) ~= Particles(i).brightestZ+1
                            Particles(i).discardThis = 1;
                            Particles(i).noIntensityAnalysis = 1;
                            falsePositives = falsePositives + 1;
                        end
                    elseif Particles(i).z(max_index -1) ~= Particles(i).brightestZ-1 ...
                        && Particles(i).z(max_index +1) ~= Particles(i).brightestZ +1
                        Particles(i).discardThis = 1;
                        Particles(i).noIntensityAnalysis = 1;
                        falsePositives = falsePositives + 1;
                    end
                elseif num_shadows == 2
                    if  Particles(i).brightestZ == Particles(i).z(end) ||...
                        Particles(i).brightestZ == Particles(i).z(1)                                
                        Particles(i).discardThis = 1;
                        Particles(i).noIntensityAnalysis = 1;
                        falsePositives = falsePositives + 1;
                    elseif Particles(i).z(max_index -1) ~= Particles(i).brightestZ-1 ...
                        || Particles(i).z(max_index +1) ~= Particles(i).brightestZ +1
                        Particles(i).discardThis = 1;
                        Particles(i).noIntensityAnalysis = 1;
                        falsePositives = falsePositives + 1;                  
                    end
                end
            end
        end

        %Create a final Spots structure to be fed into TrackmRNADynamics
        Spots{q} = [];
        fields = fieldnames(Particles);
        num_fields = length(fields);
        for i = initial_frame:num_frames
            frames = find([Particles.frame]==i);
            if ~isempty(frames)
                for j = frames(1):frames(end)
                    if ~Particles(j).discardThis
                        Spots{q}(i).Fits(j-frames(1)+1) = Particles(j);
                    end
                    %Sometimes, all Spots are discarded in a frame. In that
                    %case, create an empty Spots entry in that frame.
                    if length(Spots{q})<i
                        for l = 1:num_fields
                            Spots{q}(i).Fits.(fields{l}) = [];
                        end
                    end
                end
            else 
                for l = 1:num_fields
                    Spots{q}(i).Fits.(fields{l}) = [];
                end
            end
        end

        %Clean up Spots to remove empty rows
        Dots{q} = struct('Fits', []);
        for i = 1:length(Spots{q})
            Dots{q}(i).Fits = [];
            for j = 1:length(Spots{q}(i).Fits)
                if j~=1
                    if ~isempty(Spots{q}(i).Fits(j).z)...
                            && ~isequal(Spots{q}(i).Fits(j).CentralIntensity,Spots{q}(i).Fits(j-1).CentralIntensity)
                            Dots{q}(i).Fits = [Dots{q}(i).Fits, Spots{q}(i).Fits(j)];                    
                    end
                else
                    if ~isempty(Spots{q}(i).Fits(j).z)                        
                        Dots{q}(i).Fits = [Dots{q}(i).Fits, Spots{q}(i).Fits(j)];
                    end
                end
            end
        end
        for i = 1:length(Dots{q})
            if isstruct(Dots{q}(i).Fits)
                Spots{q}(i).Fits = rmfield(Dots{q}(i).Fits, 'r');
                Spots{q}(i).Fits = rmfield(Spots{q}(i).Fits, 'discardThis');
            else
                Spots{q}(i).Fits = [];
            end
        end

        %AR 7/10/16: Optional time tracking using track_Spots script. Also
        %makes some potentially useful plots. This was originally here to have
        %a single, fully integrated script before this segmentation was worked
        %into the rest of the pipeline.
        neighborhood = 3000 / pixelSize;
        if TrackSpots
            Particles = track_Spots(Particles, neighborhood);
            save([DropboxFolder,filesep,Prefix,filesep,'Particles_SS.mat'], 'Particles', '-v7.3');
        end
        
        %If we only have one channel, then convert Spots{q} to a
        %standard structure.
        if nCh==1
           Spots=Spots{1};
        end


        mkdir([DropboxFolder,filesep,Prefix]);
        save([DropboxFolder,filesep,Prefix,filesep,'Spots.mat'], 'Spots', '-v7.3');    
    end

    t = toc;
   disp ['Elapsed time: ',num2str(t/60),' min']
    if ~just_tifs
        logFile = [DropboxFolder,filesep,Prefix,filesep,'log.mat', '-v7.3'];
        if exist(logFile, 'file')
            load(logFile);
        else
            log = struct();
        end
        log(end+1).Date = date;
        log(end).runTime = t/60; %min
        log(end).InitialFrame = initial_frame;
        log(end).LastFrame = num_frames;
        log(end).NFrames = num_frames - initial_frame + 1;
        log(end).TimePerFrame = (t/60)/(num_frames-initial_frame + 1);

        if ~just_dog
            detectedCircles = 0;
            detectedBalls = 0;
            if iscell(Spots)
                for i = 1:length(Spots{q})
                    for j = 1:length(Spots{q}(i).Fits)
                        detectedCircles = detectedCircles + length(Spots{q}(i).Fits(j).z);
                        detectedBalls = detectedBalls + 1;
                    end
                end
            else
               for i = 1:length(Spots)
                    for j = 1:length(Spots(i).Fits)
                        detectedCircles = detectedCircles + length(Spots(i).Fits(j).z);
                        detectedBalls = detectedBalls + 1;
                    end
                end 
            end
            display(['Detected Spots: ',num2str(detectedCircles)])
            log(end).falsePositives = falsePositives;
            log(end).totalCircles = detectedCircles;
            log(end).totalBalls = detectedBalls;
            log(end).avgZSize = detectedCircles/detectedBalls;
            log(end).Threshold = Threshold;
            if isfield(log, 'Classifier')
                log(end).Classifier = log(end-1).Classifier;
            end
        else     
            log(end).Classifier = classifierPathCh1;     
        end
        save(logFile, 'log');
    end
    try
        poolobj = gcp('nocreate');
        delete(poolobj);
    catch
        %fails if the parallel pool has timed out.
    end
end
end