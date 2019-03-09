function fittedLineEquations = fitSingleTraces(Prefix,Particles,Spots,schnitzcells,FrameInfo,ElapsedTime,varargin)
% fitSingleTraces(varargin)
%
% DESCRIPTION
%
% ARGUMENTS
%
% OPTIONS
%
% Author (contact): Emma Luu (emma_luu@berkeley.edu)
% Created: 10/25/18
% Last Updated: 11/20/18
%
% Documented by: Emma Luu (emma_luu@berkeley.edu)

%% Initializing options
initialOnly = 0;
skipSavingTraces = 0;

%% Checking varargin
if ~isempty(varargin)
    for i=1:length(varargin)
        if strcmpi(varargin{i},'initalOnly')
            initialOnly = 1;
        elseif strcmpi(varargin{i},'skipSavingTraces')
            skipSavingTraces = 1;
        end
    end
end

%% Getting relevant information
%Information about about folders
[~,~,DefaultDropboxFolder,~,~]=...
    DetermineLocalFolders;

%What type of experiment are we dealing with? Get this out of MovieDatabase
[~,~,DropboxFolder,~, ~,...
    ~, ~, ~, ~, ~,~] = readMovieDatabase(Prefix);

% refactor in progress, we should replace readMovieDatabase with getExperimentDataFromMovieDatabase
[~, ~, ~, ~, ~, ~,...
    ~, ~, ~, ~, ~, ~, ~,...
    nc9, nc10, nc11, nc12, nc13, nc14, ~] = getExperimentDataFromMovieDatabase(Prefix, DefaultDropboxFolder);

%Create the particle array. This is done so that we can support multiple
%channels. Also figure out the number of channels
if iscell(Particles)
    NChannels=length(Particles);
else
    Particles={Particles};
end

if ~iscell(Spots)
    Spots={Spots};
end

%Do we need to convert any NaN chars into doubles?
if strcmpi(nc14,'nan')
    nc14=nan;
end
if strcmpi(nc13,'nan')
    nc13=nan;
end
if strcmpi(nc12,'nan')
    nc12=nan;
end
if strcmpi(nc11,'nan')
    nc11=nan;
end
if strcmpi(nc10,'nan')
    nc10=nan;
end
if strcmpi(nc9,'nan')
    nc9=nan;
end


%% Fitting shapes to single traces (includes time on and initial rate of loading)
% This section of code is will fit lines in a piece wise fashion to the
% single traces. fittedLineEquations correspond to the stored fitted lines
% of the particles, where the indexing is as follows:
% fittedLineEquations(particleNumber) with fields: Coefficients,
% ErrorEstimation, and frameIndex, which are described below. This currently
% does not support more than one channel. Please contact Emma to work on
% implementing it for two channels.

if NChannels > 1
    % implement crude? initial slope and time on calculation
else
    numberOfParticles = size(Particles{:},2);
    correspondingNCInfo = [FrameInfo.nc];
    currentChannel = 1;
    
    nuclearCycleBoundaries = [nc9,nc10,nc11,nc12,nc13,nc14]; % in units of frames
    for i = 1:length(nuclearCycleBoundaries)
        if nuclearCycleBoundaries(i) > 0
            nuclearCycleBoundaries(i) = ElapsedTime(nuclearCycleBoundaries(i)); % in units of minutes
        end
    end
    
    % The following for loop will create fittedLineEquations, which is a
    % structure array with the following fields:
    % Coefficients: a cell array of the coefficients of the lines fitted to
    % the trace, where each row is for a line. The slope has units of
    % AU/minutes and the constant has units of AU. The function polyfit was used to make this.
    % ErrorEstimation: The error estimation structure given by polyfit
    % FrameIndex: a cell array of the first and last frames used to create
    % the corredponding lines, where each row corresponds to a line.
    % Note that the first fitted line can be used to find the time on
    % (accomplished by taking the root of the equation with the function
    % roots).
    
    % Example:
    % To plot the line of the first three fitted line would be written as:
    % figure()
    % hold on
    % currentParticle = 1; % doing this for the first particle
    % for currentLine = 1:3
    %     currentCoefficients = fittedLineEquations(currentParticle).coefficients{currentLine};
    %     currentXPoints = fittedLineEquations(currentParticle).frameIndex{currentLine};
    %     currentYPoints = polyval(currentCoefficients,currentXPoints);
    %     plot(currentXPoints, currentYPoints);
    % end
    if ~skipSavingTraces
        % The folder that will save all the figures created by this function
        savedImageFolder = [DropboxFolder,filesep,Prefix,filesep,...
            'FittedInitialRates'];
        if ~exist(savedImageFolder, 'dir')
            mkdir(savedImageFolder)
        else
            delete([savedImageFolder,filesep,'*.*'])
        end
    end
    
    % lines are only fitted to traces that are longer than 3 data points
    averagingLength = 3; % average over 3 points. Current default setting
    waitBarHandle = waitbar(0,'Fitting traces');
    for currentParticle = 1:numberOfParticles
        waitbar(currentParticle/numberOfParticles,waitBarHandle);
        figure('visible','off')
        
        [frameIndex,coefficients,errorEstimation,numOfPartUsedForFit] = ...
            fitASingleTrace(currentParticle,Particles,Spots,currentChannel,...
            schnitzcells,ElapsedTime,nuclearCycleBoundaries,correspondingNCInfo,...
            averagingLength);
        
        fittedLineEquations(currentParticle).frameIndex = frameIndex;
        
        fittedLineEquations(currentParticle).Coefficients = ...
            coefficients;
        
        fittedLineEquations(currentParticle).ErrorEstimation = ...
            errorEstimation;
        
        fittedLineEquations(currentParticle).numberOfParticlesUsedForFit = ...
            numOfPartUsedForFit; % only happens for 
        
        if ~skipSavingTraces
            set(gcf,'Position',[85 43.5000 787 601]) % resizing the figure
            fileName = [savedImageFolder,filesep,...
                'Particle' num2str(currentParticle) '.png'];
            saveas(gcf,fileName)
        end
        close(gcf)
    end
    close(waitBarHandle)
end

end