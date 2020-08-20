function fillFrameGaps(prefix, varargin)
%% Information about the script
% fillFrameGaps will add the particle in frames that are
% between frames where the particle is present.

% It is strongly encouraged that you disconnect any particles that are
% actually 2 particles before running this script.

%Author: Emma Luu (emma_luu@berkeley.edu)
%% Loading the data set of interest

intScale = 1; % for intScale option
notAllNC = 0; % for desiredNC(s) option

for i = 1:length(varargin)
    if strcmp(varargin{i},'intScale')
            intScale = varargin{i+1};
    elseif strcmp(varargin{i}, 'desiredNC(s)')
            notAllNC = 1;
            desiredNC = varargin{i+1};
    end
end

display(['integration scaling factor: ', num2str(intScale)]);

% getting folder information
[~,~,dropboxFolder,~,~]= DetermineLocalFolders(prefix);
[~,~,defaultDropboxFolder,~,PreProcPath]=...
    DetermineLocalFolders;
dataFolder=[dropboxFolder,filesep,prefix];
FilePrefix=[dataFolder(length(dropboxFolder)+2:end),'_'];

% Particle Information and Loading
particlePathName = [dataFolder,filesep,'Particles.mat'];
load(particlePathName)
if ~iscell(Particles)
    Particles = {Particles};
    SpotFilter = {SpotFilter};
end
numberOfParticles = size(Particles{:},2);

NChannels=1;
nameSuffix=['_ch',iIndex(1,2)];
currentChannel = 1;

% Spots Information and Loading
spotsPathName = [dataFolder,filesep,'Spots.mat'];
load(spotsPathName)
if ~iscell(Spots)
    Spots = {Spots};
end

% FrameInfo Information and Loading
frameInfoPathName = [dataFolder,filesep,'FrameInfo.mat'];
load(frameInfoPathName)
numberOfFrames = length(FrameInfo);
[~,ExperimentType,~,~,~,~,Channel1,Channel2,~,~,~,~,~,...
    nc9, nc10,nc11,nc12,nc13,nc14,~] = getExperimentDataFromMovieDatabase(prefix,defaultDropboxFolder);
nuclearCycleBoundaries = [nc9,nc10,nc11,nc12,nc13,nc14];
nuclearCycleNumber = 9:14;
pixelSize = FrameInfo(1).PixelSize*1000; %nm
snippet_size = 2*(floor(1300/(2*pixelSize))) + 1; % nm. note that this is forced to be odd
PixelsPerLine = FrameInfo(1).PixelsPerLine;
LinesPerFrame = FrameInfo(1).LinesPerFrame;
ZSlices=FrameInfo(1).NumberSlices+2;

% Ellipses Information and Loading
try
    load([dropboxFolder,filesep,FilePrefix(1:end-1),filesep,'Ellipses.mat'])

    % Schnitzcells Information and Loading
    if ~exist([dropboxFolder,filesep,FilePrefix(1:end-1),filesep,FilePrefix(1:end-1),'_lin.mat'], 'file')
        error('Need this file!')
    end
    schnitzcellsPathName = [dropboxFolder,filesep,FilePrefix(1:end-1),filesep,FilePrefix(1:end-1),'_lin.mat'];
    load(schnitzcellsPathName)
end

% log Information and Loading
logFile = [dropboxFolder, filesep, prefix, filesep, 'log.mat'];
% assuming that the log file has already been made by segmentspots
load(logFile)

%Taken from CheckParticleTracking
%Remove the schnitz fields that can give us problems potentially if
%present. I don't know how this came to be, but it's for fields that
%are not all that relevant. The fields are: approved, ang
try
    if isfield(schnitzcells,'approved')
        schnitzcells=rmfield(schnitzcells,'approved');
    end
    if isfield(schnitzcells,'ang')
        schnitzcells=rmfield(schnitzcells,'ang');
    end
end

%% Sorting Particles by their first frames
for ChN=1:NChannels
    for currentParticle=1:length(Particles{ChN})
        FirstFrame(currentParticle)=Particles{ChN}(currentParticle).Frame(1);
    end
    [~,Permutations]=sort(FirstFrame);
    Particles{ChN}=Particles{ChN}(Permutations);
    clear FirstFrame
end

% Creating NDigits
if numberOfFrames<1E3
    NDigits=3;
elseif numberOfFrames<1E4
    NDigits=4;
end


%% Section for option desiredNC(s) 

if notAllNC
    allParticleNC = zeros(1,numberOfParticles);
    for currentParticle = 1:numberOfParticles
        [frame,~,~,~,~,~,~,~,~,~,~,~,~]=...
            GetParticleTrace(currentParticle,...
            Particles{currentChannel},Spots{currentChannel});
        correspondingNCInfo = [FrameInfo.nc];
        currentNCRange = unique(correspondingNCInfo(frame));
        allParticleNC(currentParticle) = currentNCRange(1);
    end
    
    particlesOfInterest = [];
    for currentNC = desiredNC
        particleSubset = find(allParticleNC == currentNC);
        particlesOfInterest = [particlesOfInterest particleSubset];
    end
else
    particlesOfInterest = 1:numberOfParticles; % look at all of the particles
end


%% Beginning of double checking

% Checking for particles with frames missing somewhere in the middle of
% their trace
particlesToDoubleCheck = [];
for currentParticle = particlesOfInterest
    currentFrames = Particles{currentChannel}(currentParticle).Frame;
    framesWithNoGaps = currentFrames(1):currentFrames(end);
    if ~isequal(currentFrames,framesWithNoGaps)
        particlesToDoubleCheck = [particlesToDoubleCheck currentParticle];
    end
end
disp(['Need to double check ' num2str(length(particlesToDoubleCheck)) ' particle(s).'])
disp(['Checking particles : ' num2str(particlesToDoubleCheck)])

%% Adding missing frames
% Process flow:
% 1. Go to frame where the particle is missing
% 2. Go to the previous position of the particle
% 3. Find the spots and shadows near the previous position.
% 4. Repeat until all missing frames have been filled
framesModified = {};
counter = 1;
frameInformationStruct = FrameInfo;
for currentParticle = particlesToDoubleCheck
    currentFrames = Particles{currentChannel}(currentParticle).Frame;
    framesWithNoGaps = currentFrames(1):currentFrames(end);
    framesToCheck = framesWithNoGaps(~ismember(framesWithNoGaps,currentFrames));
    framesModified{counter} = framesToCheck;
    counter = counter + 1;
    
    for currentFrame = framesToCheck
        % Need to redefine the below each time you check for a new one
        numberOfParticles = size(Particles{:},2);
        currentFrames = Particles{currentChannel}(currentParticle).Frame;
        
%         % getting current nucleus and schnitz
%         schnitzIndex=find(schnitzcells(Particles{currentChannel}(i).Nucleus).frames==currentFrame);
%         nucleusIndex=schnitzcells(Particles{currentChannel}(i).Nucleus).cellno(schnitzIndex);
%         nucleusCenterCoordinates = [Ellipses{currentFrame}(nucleusIndex,1:2)];
        
        % Finding the position in the previous frame
        previousFrame = currentFrame - 1;
        [x,y,~]=SpotsXYZ(Spots{currentChannel}(previousFrame)); % Looking at the previous frame
        currentParticleIndex=...
            Particles{currentChannel}(currentParticle).Index(currentFrames==previousFrame);
        previousPositionOfParticle = [x(currentParticleIndex) y(currentParticleIndex)];
        
        
        % The code below is from CheckParticleTracking (slightly modified)
        
        connectPositionX = previousPositionOfParticle(1);%nucleusCenterCoordinates(1); % Why do we need to add 1?
        connectPositionY = previousPositionOfParticle(2);%nucleusCenterCoordinates(2);
        
        % check that the clicked particle isn't too close to the
        % edge of the frame
        if (connectPositionX > snippet_size/2) && (connectPositionX + snippet_size/2 < PixelsPerLine)...
                && (connectPositionY > snippet_size/2) && (connectPositionY + snippet_size/2 < LinesPerFrame)
            SpotsIndex = length(Spots{currentChannel}(currentFrame).Fits)+1;
            breakflag = 0;
            parfor currentZSlice = 1:ZSlices
                spotsIm=imread([PreProcPath,filesep,FilePrefix(1:length(FilePrefix)-1),filesep,...
                    FilePrefix,iIndex(currentFrame,NDigits),'_z',iIndex(currentZSlice,2),nameSuffix,'.tif']);
                try
                    imAbove = double(imread([PreProcPath,filesep,FilePrefix(1:length(FilePrefix)-1),filesep,...
                        FilePrefix,iIndex(currentFrame,NDigits),'_z',iIndex(currentZSlice-1,2),nameSuffix,'.tif']));
                    imBelow = double(imread([PreProcPath,filesep,FilePrefix(1:length(FilePrefix)-1),filesep,...
                        FilePrefix,iIndex(currentFrame,NDigits),'_z',iIndex(currentZSlice+1,2),nameSuffix,'.tif']));
                catch
                    imAbove = nan(size(spotsIm,1),size(spotsIm,2));
                    imBelow = nan(size(spotsIm,1),size(spotsIm,2));
                end
                Threshold = min(min(spotsIm));
                dog = spotsIm;
                im_thresh = dog >= Threshold;
                [im_label, ~] = bwlabel(im_thresh);
                microscope = frameInformationStruct(1).FileMode;
                show_status = 0;
                fig = [];
                k = 1; %This is supposed to be the index for the partiles in an image.
                %However, this image only contains one particle
                neighborhood = round(1300 / pixelSize); %nm
                %Get the information about the spot on this z-slice
                
                
                temp_particles{currentZSlice} = identifySingleSpot(k, {spotsIm,imAbove,imBelow}, im_label, dog, neighborhood, snippet_size, ...
                    pixelSize, show_status, fig, microscope, [1, connectPositionX, connectPositionY], [], '', intScale);
            end
            
            % future plan: use saveParticleInformation.m
            for currentZSlice = 1:ZSlices
                if ~isempty(temp_particles{currentZSlice})
                    %Copy the information stored on temp_particles into the
                    %Spots structure
                    if ~isempty(temp_particles{currentZSlice}{1})
                        Spots{currentChannel}(currentFrame).Fits(SpotsIndex).FixedAreaIntensity(currentZSlice)=...
                            temp_particles{currentZSlice}{1}{1};
                        Spots{currentChannel}(currentFrame).Fits(SpotsIndex).xFit(currentZSlice)=...
                            temp_particles{currentZSlice}{1}{2};
                        Spots{currentChannel}(currentFrame).Fits(SpotsIndex).yFit(currentZSlice)=...
                            temp_particles{currentZSlice}{1}{3};
                        Spots{currentChannel}(currentFrame).Fits(SpotsIndex).Offset(currentZSlice)=...
                            temp_particles{currentZSlice}{1}{4};
                        Spots{currentChannel}(currentFrame).Fits(SpotsIndex).Area(currentZSlice)=...
                            temp_particles{currentZSlice}{1}{6};
                        Spots{currentChannel}(currentFrame).Fits(SpotsIndex).yDoG(currentZSlice)=...
                            temp_particles{currentZSlice}{1}{9};
                        Spots{currentChannel}(currentFrame).Fits(SpotsIndex).xDoG(currentZSlice)=...
                            temp_particles{currentZSlice}{1}{10};
                        Spots{currentChannel}(currentFrame).Fits(SpotsIndex).GaussianIntensity(currentZSlice)=...
                            temp_particles{currentZSlice}{1}{11};
                        Spots{currentChannel}(currentFrame).Fits(SpotsIndex).CentralIntensity(currentZSlice)=...
                            temp_particles{currentZSlice}{1}{12};
                        Spots{currentChannel}(currentFrame).Fits(SpotsIndex).DOGIntensity(currentZSlice)=...
                            temp_particles{currentZSlice}{1}{13};
                        Spots{currentChannel}(currentFrame).Fits(SpotsIndex).ConfidenceIntervals{currentZSlice}=...
                            temp_particles{currentZSlice}{1}{19};

                        Spots{currentChannel}(currentFrame).Fits(SpotsIndex).z(currentZSlice)=...
                            currentZSlice;
                        Spots{currentChannel}(currentFrame).Fits(SpotsIndex).gaussParams{currentZSlice}=...
                            temp_particles{currentZSlice}{1}{22};
                        Spots{currentChannel}(currentFrame).Fits(SpotsIndex).intArea=...
                            temp_particles{currentZSlice}{1}{23};
                        Spots{currentChannel}(currentFrame).Fits(SpotsIndex).discardThis=...
                            0;
                        Spots{currentChannel}(currentFrame).Fits(SpotsIndex).frame=...
                            currentFrame;
                        Spots{currentChannel}(currentFrame).Fits(SpotsIndex).r=...
                            0;
                        Spots{currentChannel}(currentFrame).Fits(SpotsIndex).FixedAreaIntensity3 = NaN;
                        Spots{currentChannel}(currentFrame).Fits(SpotsIndex).FixedAreaIntensity5 = NaN;
                        Spots{currentChannel}(currentFrame).Fits(SpotsIndex).cylIntensity(currentZSlice) = temp_particles{currentZSlice}{1}{24};
                        Spots{currentChannel}(currentFrame).Fits(SpotsIndex).brightestZ = NaN;
                        
                    else
                        Spots{currentChannel}(currentFrame).Fits(SpotsIndex).FixedAreaIntensity(currentZSlice)=...
                            nan;
                        Spots{currentChannel}(currentFrame).Fits(SpotsIndex).xFit(currentZSlice)=...
                            nan;
                        Spots{currentChannel}(currentFrame).Fits(SpotsIndex).yFit(currentZSlice)=...
                            nan;
                        Spots{currentChannel}(currentFrame).Fits(SpotsIndex).Offset(currentZSlice)=nan;  
                        Spots{currentChannel}(currentFrame).Fits(SpotsIndex).Are(currentZSlice)...
                            nan;
                        Spots{currentChannel}(currentFrame).Fits(SpotsIndex).yDoG(currentZSlice)=...
                            nan;
                        Spots{currentChannel}(currentFrame).Fits(SpotsIndex).xDoG(currentZSlice)=...
                            nan;
                        Spots{currentChannel}(currentFrame).Fits(SpotsIndex).GaussianIntensity(currentZSlice)=...
                            nan;
                        Spots{currentChannel}(currentFrame).Fits(SpotsIndex).CentralIntensity(currentZSlice)=...
                            nan;
                        Spots{currentChannel}(currentFrame).Fits(SpotsIndex).DOGIntensity(currentZSlice)=...
                            nan;

                        Spots{currentChannel}(currentFrame).Fits(SpotsIndex).ConfidenceIntervals{currentZSlice}=...
                            nan;
                        Spots{currentChannel}(currentFrame).Fits(SpotsIndex).z(currentZSlice)=...
                            currentZSlice;
                        Spots{currentChannel}(currentFrame).Fits(SpotsIndex).gaussParams{currentZSlice}=...
                            nan;
                        Spots{currentChannel}(currentFrame).Fits(SpotsIndex).discardThis=...
                            0;
                        Spots{currentChannel}(currentFrame).Fits(SpotsIndex).frame=...
                            currentFrame;
                        Spots{currentChannel}(currentFrame).Fits(SpotsIndex).r=...
                            0;
                        Spots{currentChannel}(currentFrame).Fits(SpotsIndex).FixedAreaIntensity3 = NaN;
                        Spots{currentChannel}(currentFrame).Fits(SpotsIndex).cylIntensity(currentZSlice) = NaN;
                        Spots{currentChannel}(currentFrame).Fits(SpotsIndex).FixedAreaIntensity5 = NaN;
                        Spots{currentChannel}(currentFrame).Fits(SpotsIndex).brightestZ = NaN;
                    end
                else
                    disp('No spot added. Did you click too close to the image boundary?')
                    breakflag = 1;
                    break
                end
            end
            
            if ~breakflag
                force_z = 0;
                use_integral_center = 1;
                [tempSpots,~] = findBrightestZ(Spots{currentChannel}(currentFrame).Fits(SpotsIndex), -1, use_integral_center, force_z);
                Spots{currentChannel}(currentFrame).Fits(SpotsIndex) = tempSpots;
                
                %Add this to SpotFilter, which tells the code that this spot is
                %above the threshold. First, check whether the
                %dimensions of SpotFilter need to be altered. If so, pad it with NaNs
                if size(SpotFilter{currentChannel},2)>SpotsIndex
                    SpotFilter{currentChannel}(currentFrame,SpotsIndex)=1;
                else
                    %Pad with NaNs
                    SpotFilter{currentChannel}(:,end:SpotsIndex)=NaN;
                    SpotFilter{currentChannel}(currentFrame,SpotsIndex)=1;
                end
                %
                %Turn this spot into a new particle. This is the equivalent of
                %the 'u' command.
                [SpotFilter{currentChannel},Particles{currentChannel}]=...
                    TransferParticle(Spots{currentChannel},...
                    SpotFilter{currentChannel},Particles{currentChannel},...
                    currentFrame,SpotsIndex);
                numberOfParticles = numberOfParticles + 1;
                Particles{currentChannel}(end).FrameApproved = true;
                
                %Connect this particle to the CurrentParticle. This is
                %the equivalent of running the 'c' command.
                
                Particles{currentChannel}=...
                    JoinParticleTraces(currentParticle,...
                    numberOfParticles,Particles{currentChannel});
                
                %Finally, force the code to recalculate the fluorescence trace
                %for this particle
                PreviousParticle=0;
                spotAddedMessage = ['Spot added to the current particle ' num2str(currentParticle)];
                disp(spotAddedMessage)
            else
                warning('This particle is too close to the edge. A spot can''t be added here.');
            end
        end
    end
end

%% Saving Particles and Spots
% Possible insert: Double check the particles' traces before asking to save them
log(end+1).Date = date;
log(end).integrationScale = intScale;
log(end).FunctionOrScriptUsed = 'addMissingFramesToParticles';
log(end).ParticlesEdited = particlesToDoubleCheck;
log(end).CorrespondingFramesAdded = framesModified;

% usersAnswer = input('Would you like to save these traces? (y/n)','s');
% saving = isequal('y',lower(usersAnswer));
saving = true;
if saving && ~isempty(particlesToDoubleCheck)
    disp('Saving the changes')
    % Saving all that hardwork!

    %If we only have one channel bring Particles back to the legacy
    %format without any cells

    if NChannels==1
        Particles=Particles{1};
        Spots=Spots{1};
        SpotFilter=SpotFilter{1};
    end
     
    save(particlePathName,'Particles','SpotFilter', '-v7.3')
    save(spotsPathName,'Spots','-v7.3')
    save(logFile,'log','-v7.3')
    % Note: did not use histone overlay so the _lin.mat file was not saved.
end

end