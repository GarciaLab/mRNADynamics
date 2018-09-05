function fillGaps(prefix, varargin)
%% Information about the script
% fillGaps will add the particle in frames that are
% between frames where the particle is present.

% It is strongly encouraged that you disconnect any particles that are
% actually 2 particles before running this script.

%Author: Emma Luu (emma_luu@berkeley.edu)
%% Loading the data set of interest

%[prefix,~] = getPrefixAndFolder; % in case the prefix above is not yours.
intScale = 1;

for i = 1:length(varargin)
    if strcmpi(varargin{i}, 'intScale')
        intScale = varargin{i+1};
        display(['integration scaling factor: ', num2str(intScale)]);
    end
end

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
    for i=1:length(Particles{ChN})
        FirstFrame(i)=Particles{ChN}(i).Frame(1);
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
else
    error('No more than 10,000 frames supported. Change this in the code')
end

%% Beginning of double checking

% Checking for particles with frames missing somewhere in the middle of
% their trace
particlesToDoubleCheck = [];
for i = 1:numberOfParticles
    currentFrames = Particles{currentChannel}(i).Frame;
    framesWithNoGaps = currentFrames(1):currentFrames(end);
    if ~isequal(currentFrames,framesWithNoGaps)
        particlesToDoubleCheck = [particlesToDoubleCheck i];
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
for i = particlesToDoubleCheck
    currentFrames = Particles{currentChannel}(i).Frame;
    framesWithNoGaps = currentFrames(1):currentFrames(end);
    framesToCheck = framesWithNoGaps(~ismember(framesWithNoGaps,currentFrames));
    framesModified{counter} = framesToCheck;
    counter = counter + 1;
    
    for currentFrame = framesToCheck
        % Need to redefine the below each time you check for a new one
        numberOfParticles = size(Particles{:},2);
        currentFrames = Particles{currentChannel}(i).Frame;
        
        % getting current nucleus and schnitz
        schnitzIndex=find(schnitzcells(Particles{currentChannel}(i).Nucleus).frames==currentFrame);
        nucleusIndex=schnitzcells(Particles{currentChannel}(i).Nucleus).cellno(schnitzIndex);
        nucleusCenterCoordinates = [Ellipses{currentFrame}(nucleusIndex,1:2)];
        
        % Finding the position in the previous frame
        previousFrame = currentFrame - 1;
        [x,y,z]=SpotsXYZ(Spots{currentChannel}(previousFrame)); % Looking at the previous frame
        currentParticleIndex=...
            Particles{currentChannel}(i).Index(currentFrames==previousFrame);
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
            parfor j = 1:ZSlices
                spotsIm=imread([PreProcPath,filesep,FilePrefix(1:length(FilePrefix)-1),filesep,...
                    FilePrefix,iIndex(currentFrame,NDigits),'_z',iIndex(j,2),nameSuffix,'.tif']);
                try
                    imAbove = double(imread([PreProcPath,filesep,FilePrefix(1:length(FilePrefix)-1),filesep,...
                        FilePrefix,iIndex(currentFrame,NDigits),'_z',iIndex(j-1,2),nameSuffix,'.tif']));
                    imBelow = double(imread([PreProcPath,filesep,FilePrefix(1:length(FilePrefix)-1),filesep,...
                        FilePrefix,iIndex(currentFrame,NDigits),'_z',iIndex(j+1,2),nameSuffix,'.tif']));
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
                
                
                temp_particles{j} = identifySingleSpot(k, {spotsIm,imAbove,imBelow}, im_label, dog, neighborhood, snippet_size, ...
                    pixelSize, show_status, fig, microscope, [1, connectPositionX, connectPositionY], [], '', intScale);
            end
            
            % future plan: use saveParticleInformation.m
            for j = 1:ZSlices
                if ~isempty(temp_particles{j})
                    %Copy the information stored on temp_particles into the
                    %Spots structure
                    if ~isempty(temp_particles{j}{1})
                        Spots{currentChannel}(currentFrame).Fits(SpotsIndex).FixedAreaIntensity(j)=...
                            temp_particles{j}{1}{1};
                        Spots{currentChannel}(currentFrame).Fits(SpotsIndex).xFit(j)=...
                            temp_particles{j}{1}{2};
                        Spots{currentChannel}(currentFrame).Fits(SpotsIndex).yFit(j)=...
                            temp_particles{j}{1}{3};
                        Spots{currentChannel}(currentFrame).Fits(SpotsIndex).Offset(j)=...
                            temp_particles{j}{1}{4};
                        Spots{currentChannel}(currentFrame).Fits(SpotsIndex).Area{j}=...
                            temp_particles{j}{1}{6};
                        Spots{currentChannel}(currentFrame).Fits(SpotsIndex).yDoG(j)=...
                            temp_particles{j}{1}{9};
                        Spots{currentChannel}(currentFrame).Fits(SpotsIndex).xDoG(j)=...
                            temp_particles{j}{1}{10};
                        Spots{currentChannel}(currentFrame).Fits(SpotsIndex).GaussianIntensity(j)=...
                            temp_particles{j}{1}{11};
                        Spots{currentChannel}(currentFrame).Fits(SpotsIndex).CentralIntensity(j)=...
                            temp_particles{j}{1}{12};
                        Spots{currentChannel}(currentFrame).Fits(SpotsIndex).DOGIntensity(j)=...
                            temp_particles{j}{1}{13};
                        Spots{currentChannel}(currentFrame).Fits(SpotsIndex).ConfidenceIntervals{j}=...
                            temp_particles{j}{1}{19};

                        Spots{currentChannel}(currentFrame).Fits(SpotsIndex).z(j)=...
                            j;
                        Spots{currentChannel}(currentFrame).Fits(SpotsIndex).gaussParams{j}=...
                            temp_particles{j}{1}{22};
                        Spots{currentChannel}(currentFrame).Fits(SpotsIndex).intArea=...
                            temp_particles{j}{1}{23};
                        Spots{currentChannel}(currentFrame).Fits(SpotsIndex).discardThis=...
                            0;
                        Spots{currentChannel}(currentFrame).Fits(SpotsIndex).frame=...
                            currentFrame;
                        Spots{currentChannel}(currentFrame).Fits(SpotsIndex).r=...
                            0;
                        Spots{currentChannel}(currentFrame).Fits(SpotsIndex).FixedAreaIntensity3 = NaN;
                        Spots{currentChannel}(currentFrame).Fits(SpotsIndex).FixedAreaIntensity5 = NaN;
                        Spots{currentChannel}(currentFrame).Fits(SpotsIndex).cylIntensity(j) = temp_particles{j}{1}{24};
                        Spots{currentChannel}(currentFrame).Fits(SpotsIndex).brightestZ = NaN;
                        
                    else
                        Spots{currentChannel}(currentFrame).Fits(SpotsIndex).FixedAreaIntensity(j)=...
                            nan;
                        Spots{currentChannel}(currentFrame).Fits(SpotsIndex).xFit(j)=...
                            nan;
                        Spots{currentChannel}(currentFrame).Fits(SpotsIndex).yFit(j)=...
                            nan;
                        Spots{currentChannel}(currentFrame).Fits(SpotsIndex).Offset(j)=nan;  
                        Spots{currentChannel}(currentFrame).Fits(SpotsIndex).Area{j}=...
                            nan;
                        Spots{currentChannel}(currentFrame).Fits(SpotsIndex).yDoG(j)=...
                            nan;
                        Spots{currentChannel}(currentFrame).Fits(SpotsIndex).xDoG(j)=...
                            nan;
                        Spots{currentChannel}(currentFrame).Fits(SpotsIndex).GaussianIntensity(j)=...
                            nan;
                        Spots{currentChannel}(currentFrame).Fits(SpotsIndex).CentralIntensity(j)=...
                            nan;
                        Spots{currentChannel}(currentFrame).Fits(SpotsIndex).DOGIntensity(j)=...
                            nan;

                        Spots{currentChannel}(currentFrame).Fits(SpotsIndex).ConfidenceIntervals{j}=...
                            nan;
                        Spots{currentChannel}(currentFrame).Fits(SpotsIndex).z(j)=...
                            j;
                        Spots{currentChannel}(currentFrame).Fits(SpotsIndex).gaussParams{j}=...
                            nan;
                        Spots{currentChannel}(currentFrame).Fits(SpotsIndex).discardThis=...
                            0;
                        Spots{currentChannel}(currentFrame).Fits(SpotsIndex).frame=...
                            currentFrame;
                        Spots{currentChannel}(currentFrame).Fits(SpotsIndex).r=...
                            0;
                        Spots{currentChannel}(currentFrame).Fits(SpotsIndex).FixedAreaIntensity3 = NaN;
                        Spots{currentChannel}(currentFrame).Fits(SpotsIndex).cylIntensity(j) = NaN;
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
                    JoinParticleTraces(i,...
                    numberOfParticles,Particles{currentChannel});
                
                %Finally, force the code to recalculate the fluorescence trace
                %for this particle
                PreviousParticle=0;
                spotAddedMessage = ['Spot added to the current particle ' num2str(i)];
                disp(spotAddedMessage)
            else
                warning('This particle is too close to the edge. A spot can''t be added here.');
            end
        end
    end
end
% Need to redefine it again in case numberOfParticles is used later...
% numberOfParticles = size(Particles{:},2);

%% Saving Particles and Spots
% Possible insert: Double check the particles' traces before asking to save them
log(end+1).Date = date;
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
     
    save(particlePathName,'Particles','SpotFilter','Threshold1','Threshold2', '-v7.3')
    save(spotsPathName,'Spots','-v7.3')
    save(logFile,'log','-v7.3')
    % Note: did not use histone overlay so the _lin.mat file was not saved.
end

end