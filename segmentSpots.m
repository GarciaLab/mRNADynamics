function segmentSpots(Prefix,Threshold,varargin)
% segmentSpots(Prefix, Threshold, [Options])
%
% DESCRIPTION
% Identify and segment individual transcription spots.
%
% ARGUMENTS
% Prefix: Prefix of the data set to analyze
% Threshold: Threshold to be used. Should be kept at ~90-200 for lattice
%           light-sheet data, and at ~5-10 for confocal data (Leica SP8).
%           If left empty, then the code just generates the DoG files.
% [Options]: See below.
%
% OPTIONS
% 'displayFigures':   If you want to display plots and images.
%                
% 'TrackSpots':   Do you want to use this code to track the particles instead
%                of using TrackmRNADynamics? 
% 'Frames',N:     Run the code from frame 1 to frame N. Defaults to all
%                frames. It's suggested to run 5-20 frames for debugging.
% 'customSigmas': Prompts you to enter your custom sigmas to do DoG
%                 filtering with
% 'Shadows':    	 This option should be followed by 0, 1 or 2. This
%                specifies the number of requisite z-planes above and/or below the
%                brightest plane for a spot to have to pass quality control. 
% 'customFilters': Choose which filter to use to segment the image. Name
%                  should be a string, followed by a cell with your filter
%                  or filters
%                 ex. segmentSpots(Prefix,[],'customFilter','Structure_largest',{1,8})
%           Filter Options:
%               'Gaussian_blur'             'Median'
%               'Edges'                     'Maximum'
%               'Laplacian'                 'Minimum'
%               'Mean'                      'Std'
%               'Hessian_largest'           'Hessian_smallest'
%               [DEFAULT] 'Difference_of_Gaussian' (2 sigmas) [DEFAULT]
%               'Structure_largest' (2 sigmas)
%               'Structure_smallest' (2 sigmas)
%               
% OUTPUT
% 'Spots':  A structure array with a list of detected transcriptional loci
% in each frame and their properties.
%
% Author (contact): Armando Reimer (areimer@berkeley.edu)
% Created: 01/01/2016
% Last Updated: 12/31/2016
%
% Documented by: Armando Reimer (areimer@berkeley.edu)

warning('off','MATLAB:MKDIR:DirectoryExists');

%Default options
displayFigures=0;
trackSpots=0;
numFrames=0;
numShadows = 2;
customSigmas = 0;
customFilter = 0;
filterType = 'Difference_of_Gaussian';

for i=1:length(varargin)
    if strcmp(varargin{i},'displayFigures')
        displayFigures=1;
    elseif strcmp(varargin{i},'TrackSpots')
        trackSpots=1;
    elseif strcmp(varargin{i},'Shadows')
        if ~isnumeric(varargin{i+1}) || varargin{i+1} > 2
            error('Wrong input parameters. After ''Shadows'' you should input number of shadows (0, 1 or 2)')
        else
            numShadows=varargin{i+1};
        end
    elseif strcmpi(varargin{i}, 'customSigmas')
        customSigmas = 1;
    elseif strcmp(varargin{i},'Frames')
        if ~isnumeric(varargin{i+1})
            error('Wrong input parameters. After ''Frames'' you should input the number of frames')
        else
            numFrames=varargin{i+1};
        end
    elseif strcmp(varargin{i},'customFilter')
        customFilter = 1;
        try
            filterType = varargin{i+1};
        catch
            warning('Entered filter not recognized. Defaulting to DoG')
        end
        if iscell(varargin{i+2})
            sigmas = varargin{i+2};
            if strcmp(filterType,'Difference_of_Gaussian') || ...
                    strcmp(filterType,'Structure_largest') || ...
                    strcmp(filterType,'Structure_smallest')
                if length(sigmas) ~= 2
                    error('DoG and Structure filters require two sigma values e.g.{lower_sigma,higher_sigma}')
                end
            else
                if length(sigmas) ~= 1
                    error('All filters besides DoG and Structure require only 1 sigma value')
                end
            end
        else
            error('Entered sigma(s) not recognized. Make sure the sigma(s) are entered as numbers in a cell {}')
        end
    end
end

%If no threshold was specified, then just generate the DoG images
justDoG=0;
try
    if isempty(Threshold)
        justDoG=1;
    end
catch
    error('Please pass the argument "[]" to generate DoG images')
end
%%
tic;

maxWorkers = 56;
try
    parpool(maxWorkers);  % 6 is the number of cores the Garcia lab server can reasonably handle per user at present.
catch
    try
        parpool; %in case there aren't enough cores on the computer 
    catch
    end
    %parpool throws an error if there's a pool already running. 
end

[~,~,~,~,~,~,~,ExperimentType, Channel1, Channel2,~] =...
    readMovieDatabase(Prefix);

[SourcePath,FISHPath,DropboxFolder,MS2CodePath,PreProcPath]=...
    DetermineLocalFolders(Prefix);

load([DropboxFolder,filesep,Prefix,filesep,'FrameInfo.mat']);
microscope = FrameInfo(1).FileMode;
zSize = FrameInfo(1).NumberSlices + 2;
if numFrames == 0
    numFrames = length(FrameInfo);
end
OutputFolder1=[FISHPath,filesep,Prefix,'_',filesep,'dogs'];
mkdir(OutputFolder1) 

nCh = 1;
if strcmpi(ExperimentType, '2spot2color')
    nCh = 2;
end

%Load flat-field. We need to process this file differently the images come
%from a laser scanning or spinning disk microscope.
doFF = 1;
try
    ffim = imread([PreProcPath, filesep, Prefix, filesep,Prefix,'_FF.tif']);    
    %If we have a spinning disk confocal
    if strcmpi(FrameInfo(1).FileMode,'dspin')
        warning('Assuming a spinning disk confocal for flat-field correction')
        %Note that we brought this back to the same parameters as for a
        %LSC. We need to figure out what's going on with the flatfields on
        %the spinning disk. If not, we can always crop the image.
        ffim = CPsmooth(ffim,'Gaussian Filter',256,0);
	%If not, we assume we have a laser scanning confocal
    else
        warning('Assuming a laser-scanning confocal for flat-field correction')
        ffim = CPsmooth(ffim,'Gaussian Filter',256,0);
    end
    %Normalize the image
    ffim = double(ffim)/double(max(max(ffim)));
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

           
all_frames = cell(numFrames, zSize);
close all force;

% Support for inputoutput mode (since the coatChannel might not be
% channel1 in inputoutput ExperimentType (YJK : 1/11/2018)
%(MT, 2018-02-11) Added support for lattice imaging, maybe temporary -
%FIX LATER
if strcmpi(ExperimentType,'inputoutput') || strcmpi(ExperimentType,'lattice')
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

if justDoG
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Generate difference of Gaussian images if no threshold was given.

filterSize = round(2000/pixelSize);

    %If customSigma is desired, prompts for the sigma1 and sigma2 values
    if customSigmas
        prompt = {'Enter Sigma1 (Note: Sigma1<Sigma2):','Enter Sigma2:'};
        dlg_title = 'Input';
        num_lines = 1;
        defaultans = {'1','200'};
        answer = inputdlg(prompt,dlg_title,num_lines,defaultans);
        sigma1 = str2num(answer{1,1});
        sigma2 = str2num(answer{2,1});
    elseif ~customFilter
        %Initialize Difference of Gaussian filter parameters. filterSize >> sigma2
        %> sigma1
        sigma1 = pixelSize / pixelSize; %width of narrower Gaussian
        sigma2 = round(42000 / pixelSize); % width of wider Gaussian. AR 1/10/18: what is this number.
        sigmas = {sigma1,sigma2};
    end   
    
     
    for q = 1:nCh
%         h=waitbar(0,'Generating DoG images');
        %(MT, 2018-02-11) Added support for lattice imaging, maybe 
        %temporary - FIX LATER
        if strcmpi(ExperimentType,'inputoutput') || strcmpi(ExperimentType,'lattice')
            nameSuffix= ['_ch',iIndex(coatChannel,2)];
        else
            nameSuffix = ['_ch',iIndex(q,2)];
        end
        
        for current_frame = 1:numFrames
%             waitbar(current_frame/numFrames,h);
            if displayFigures
                for i = 1:zSize
                    im = double(imread([PreProcPath,filesep,Prefix,filesep,Prefix,'_',iIndex(current_frame,3),'_z',iIndex(i,2),nameSuffix,'.tif']));
                    if strcmp(filterType,'Difference_of_Gaussian')
                        dog = filterImage(im,filterType,sigmas, filterSize);
                    else
                        dog = filterImage(im,filterType,sigmas, []) + 100;
                    end
                    dog = padarray(dog(filterSize:end-filterSize-1, filterSize:end-filterSize-1), [filterSize,filterSize]);
                    dog_name = ['DOG_',Prefix,'_',iIndex(current_frame,3),'_z',iIndex(i,2),nameSuffix,'.tif'];
                    imwrite(uint16(dog), [OutputFolder1,filesep,dog_name])
                    imshow(dog,[]);
                end
            else 
                parfor i = 1:zSize    
                    im = double(imread([PreProcPath,filesep,Prefix,filesep,Prefix,'_',iIndex(current_frame,3),'_z',iIndex(i,2),nameSuffix,'.tif']));
                    if strcmp(filterType,'Difference_of_Gaussian')
                        dog = filterImage(im,filterType,sigmas, filterSize);
                    else
                        dog = filterImage(im,filterType,sigmas, [])+100;
                    end
                    dog = padarray(dog(filterSize:end-filterSize-1, filterSize:end-filterSize-1), [filterSize,filterSize]);
                    dog_name = ['DOG_',Prefix,'_',iIndex(current_frame,3),'_z',iIndex(i,2),nameSuffix,'.tif'];
                    imwrite(uint16(dog), [OutputFolder1,filesep,dog_name])
                end
            end
        end
%         close(h);
    end 
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Segment transcriptional loci
else       
    thresh = Threshold; %copy so we can change the value of Threshold for each channel iteration
    for q = 1:nCh
        %(MT, 2018-02-11) Added support for lattice imaging, maybe 
        %temporary - FIX LATER
        if strcmpi(ExperimentType,'inputoutput') ||  strcmpi(ExperimentType,'lattice')
            nameSuffix= ['_ch',iIndex(coatChannel,2)];
        else
            nameSuffix = ['_ch',iIndex(q,2)];
        end
        
        h=waitbar(0,'Segmenting spots');
        Threshold = thresh(q);
        for current_frame = 1:numFrames
            w = waitbar(current_frame/numFrames,h);
            set(w,'units', 'normalized', 'position',[0.4, .15, .25,.1]);
            for i = 1:zSize   
                im = double(imread([PreProcPath,filesep,Prefix,filesep,Prefix,'_',iIndex(current_frame,3),'_z',iIndex(i,2),nameSuffix,'.tif']));
                try
                    dog = double(imread([OutputFolder1, filesep,'DOG_',Prefix,'_',iIndex(current_frame,3),'_z',iIndex(i,2),nameSuffix,'.tif']));
                catch
                    error('Please re-run with threshold ''[]'' to create DoG files')
                end
                if displayFigures
                    fig = figure(1);
                    imshow(dog,[]);
                else
                    fig=[];
                end
                %apply flatfield correction
                if doFF && sum(size(im)==size(ffim))
                    im = im./ffim;
                end
                %
                im_thresh = dog >= Threshold;
                [im_label, n_spots] = bwlabel(im_thresh); 
                centroids = regionprops(im_thresh, 'centroid');


                temp_frames = {};
                temp_particles = cell(1, n_spots);
                
                if n_spots ~= 0
                    if ~displayFigures
                        parfor k = 1:n_spots
                            centroid = round(centroids(k).Centroid);
                            temp_particles(k) = identifySingleSpot(k, im, im_label, dog, ...
                                neighborhood, snippet_size, pixelSize, displayFigures, fig, microscope, 0, centroid, '');
                        end
                    else
                        for k = 1:n_spots
                            centroid = round(centroids(k).Centroid);
                            temp_particles(k) = identifySingleSpot(k, im, im_label, dog, ...
                                neighborhood, snippet_size, pixelSize, displayFigures, fig, microscope, 0, centroid, '');
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

    %%
    %Create a useful structure that can be fed into pipeline
     if ~justDoG
        Particles = struct;
        n = 1;
        h=waitbar(0,'Saving particle information');
        for i = 1:numFrames  
            waitbar(i/numFrames,h)
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
            neighborhood = round(1300 / pixelSize);
            for n = 1:numFrames
                waitbar(n/numFrames,h)
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
            if trackSpots
                for j = 1:numel(fields)-2 %do not include fields 'r' or 'frame'
                    Particles(i).(fields{j}) = Particles(i).(fields{j})(max_index);
                end
            else
            Particles(i).brightestZ = Particles(i).z(max_index);       
                if numShadows == 1
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
                elseif numShadows == 2
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
        for i = 1:numFrames
            frames = find([Particles.frame]==i);
            if ~isempty(frames)
                for j = frames(1):frames(end)
                    if ~Particles(j).discardThis
                        Spots{q}(i).Fits(j-frames(1)+1) = Particles(j);
                    end
                    %Sometimes, all spots are discarded in a frame. In that
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

        %AR 7/10/16: Optional time tracking using track_spots script. Also
        %makes some potentially useful plots. This was originally here to have
        %a single, fully integrated script before this segmentation was worked
        %into the rest of the pipeline.
        if trackSpots
            neighborhood = round(3000 / pixelSize);
            Particles = track_spots(Particles, neighborhood, numFrames);
            save([DropboxFolder,filesep,Prefix,filesep,'Particles_SS.mat'], 'Particles');
        end
        %If we only have one channel, then convert Spots to a
        %standard structure.
        if nCh==1
           Spots=Spots{1};
        end

        mkdir([DropboxFolder,filesep,Prefix]);
        save([DropboxFolder,filesep,Prefix,filesep,'Spots.mat'], 'Spots','-v7.3');    

    end

    t = toc;
    display(['Elapsed time: ',num2str(t/60),' min'])

    if ~justDoG
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
    end
end
try
    poolobj = gcp('nocreate');
    delete(poolobj);
catch
end
end