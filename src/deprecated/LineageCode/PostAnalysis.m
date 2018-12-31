
function [lrep, furrow]=PostAnalysis(Prefix,Direction,delay)
% Performs all the post analysis from an analysed data set. Run this from a
% folder where you want to save all the data. Typically in your dropbox
% folder under a different sub-folder
%% Extract all data
%Get the folders, including the default Dropbox one
[SourcePath, FISHPath, DefaultDropboxFolder, DropboxFolder, MS2CodePath, PreProcPath,...
configValues, movieDatabasePath] = DetermineAllLocalFolders(Prefix);

if ~exist('Direction','var')
    Direction=input('No Direction of Repression was indicated. Please indicate a direction now (u for up, or d for down)');
end


%Determine division times
%Load the information about the nc from moviedatabase file
[Date, ExperimentType, ExperimentAxis, CoatProtein, StemLoop, APResolution,...
Channel1, Channel2, Objective, Power, DataFolder, DropboxFolderName, Comments,...
nc9, nc10, nc11, nc12, nc13, nc14, CF] = getExperimentDataFromMovieDatabase(Prefix, DefaultDropboxFolder)

if strcmp(ExperimentAxis,'AP')
    error('This seems to be an experiment on the AP axis. Either check your prefix or change the respective MovieDatabase entry');
end


%Find the corresponding entry in moviedatabase file
if (~isempty(findstr(Prefix,'Bcd')))&(isempty(findstr(Prefix,'BcdE1')))&...
        (isempty(findstr(Prefix,'NoBcd')))&(isempty(findstr(Prefix,'Bcd1x')))&(isempty(findstr(Prefix,'Bcd4x')))
    warning('This step in CheckParticleTracking will most likely have to be modified to work')
end

% Extract Ellipses
try
    dummyEllipses=load([DropboxFolder,filesep,Prefix,filesep,'Ellipses.mat']);
    Ellipses=dummyEllipses.Ellipses;
catch
    error('No Ellipses found!');
end

frames=size(Ellipses,1);


%This checks whether all ncs have been defined
ncCheck=[nc9,nc10,nc11,nc12,nc13,nc14];
if length(ncCheck)~=6
    error('Check the nc frames in the MovieDatabase entry. Some might be missing')
end

%Do we need to convert any NaN chars into doubles?
if strcmp(lower(nc14),'nan')
    nc14=nan;
end
if strcmp(lower(nc13),'nan')
    nc13=nan;
end
if strcmp(lower(nc12),'nan')
    nc12=nan;
end
if strcmp(lower(nc11),'nan')
    nc11=nan;
end
if strcmp(lower(nc10),'nan')
    nc10=nan;
end
if strcmp(lower(nc9),'nan')
    nc9=nan;
end
if strcmp(lower(cf),'nan')
    cf=nan;
end

ncs=[nc9,nc10,nc11,nc12,nc13,nc14];

if (length(find(isnan(ncs)))==length(ncs))|(length(ncs)<6)
    error('Have the ncs been defined in MovieDatabase?')
end

%Now do the nuclear segmentation and lineage tracking. This should be put
%into an independent function.

%Create the cell array with the names.
D=dir([PreProcPath,filesep,Prefix,filesep,'*His*.tif']);
for i=1:length(D)
    names{i}=[PreProcPath,filesep,Prefix,filesep,D(i).name];
end

%% Extract images and structures
DataFolder=[DropboxFolder,filesep,Prefix];
FilePrefix=[DataFolder(length(DropboxFolder)+2:end),'_'];

% Load the lineage files
try
    dummyLineage=load([DropboxFolder,filesep,Prefix,filesep,Prefix,'_lin.mat']);
    schnitzcells=dummyLineage.schnitzcells;
catch
    error('Unable to load the lineage file. Cannot Proceed');
end

% Load the CompiledParticles
CompiledParticles=struct;
FrameInfo=struct;
cpexist=true; % Does Compiled Particles Exist?
try
    dummyParticles=load([DropboxFolder,filesep,Prefix,filesep,'CompiledParticles.mat']);
    CompiledParticles=dummyParticles.CompiledParticles;
catch
    warning('No Particles.mat found. No changes will be made to the Particles structure.');
    cpexist=false;
end
try
    dummyFrameInfo=load([DropboxFolder,filesep,Prefix,filesep,'FrameInfo.mat']);
    FrameInfo=dummyFrameInfo.FrameInfo;
catch
    error('No Frame Info found!');
end


dummyFrame=FrameInfo(1);
dimy=dummyFrame.LinesPerFrame;
dimx=dummyFrame.PixelsPerLine;

%% Initialisation
close all;
% We first get all the required variables
if Direction~='u' && Direction~='d' && Direction~='b'
    warning('The direction specified is not valid. We assume repression happens in the upward direction');
    Direction='u';   % 'u' if repression happens in the upward direction and 'd'
    % if repression happens in the downwards direction
end
if isnan(cf)
    cf=frames;
end
if nc12==0
    nc12=1;
end
mkdir([DropboxFolder,filesep,Prefix,filesep,'PostAnalyses']);
cd([DropboxFolder,filesep,Prefix,filesep,'PostAnalyses']);
%% Finding the ectogenic boundary
if cpexist
    [ecto,expline,alltracesy,alltraces,er,runs]=FindEctogenicBoundary(CompiledParticles,frames,Direction,nc12,dimy,dimx,cf);
else
    warning('No Compiled Particles found, so no Transciptional Boundary found');
end

%% Midline
 if exist('PostAnalysisWorkspace.mat')>0
     load('PostAnalysisWorkspace');
 else
    if ~isnan(cf)
        midlines=FindMidline(PreProcPath, Prefix, Ellipses, cf,dimx,ecto);
        hold off;
        figure;
        
        % Plot the distance between the morphological feature, the ventral
        % midline, and the transcriptional feature, the boundary between the
        % mesoderm and the neurogenic ectoderm
        numofframes=size(Ellipses,1);
        plot(cf:numofframes, midlines(cf:frames)-expline(cf:frames)')
        title('Distance between midline and Current Limit of Expression');
        xlabel('Frames')
        ylabel('Position along DV');
        saveas(gcf,'midline_expr_dist.png');
        
        plot(midlines(cf:frames)-ecto(2)*ones(frames-cf+1,1))
        title('Distance between Midline and Ectogenic Boundary')
        xlabel('Frames')
        ylabel('Position along DV');
        saveas(gcf,'midline_ecto_dist.png');
        hold off;
        
        diff=mean(midlines(cf:frames))-ecto(1);
    end
end

StatPlots(ecto,alltracesy,mean(midlines(cf:frames)),numofframes,nc12,alltraces,Ellipses,er)
%DummyPlots;
%% TwoSpot
save('PostAnalysisWorkspace') % save MATLAB workspace
try
    combos=twospot(Prefix,Ellipses,schnitzcells,CompiledParticles,PreProcPath, nc14,nc12,delay);
    save('PostAnalysisWorkspace') % save MATLAB workspace
 catch
%     display('Sorry Emilia, the picture didnt come out. I will fix it after I come back on Monday');
 end
