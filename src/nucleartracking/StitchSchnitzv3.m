function [schnitzcells, Ellipses] = StitchSchnitz(Prefix, nWorkers)

%This function joins schnitzcells that overlap in space and are contiguous in time.

%% load stuff
liveExperiment = LiveExperiment(Prefix);

DropboxFolder = liveExperiment.resultsFolder;
anaphaseFrames = liveExperiment.anaphaseFrames';

nc9 = anaphaseFrames(1);
nc10 = anaphaseFrames(2);
nc11 = anaphaseFrames(3);
nc12 = anaphaseFrames(4);
nc13 = anaphaseFrames(5); 
nc14 = anaphaseFrames(6);


ncVector=[0,0,0,0,0,0,0,0,nc9,nc10,nc11,nc12,nc13,nc14];

Ellipses = getEllipses(liveExperiment);
schnitzcells = getSchnitzcells(liveExperiment); 

%% setup stuff

[nFrames,~] = size(Ellipses);
RadiusPerFrame = []; %a vector of length = frames that will contain average ellipse radii per frame:
for fr = 1:nFrames
    if ~isempty(Ellipses{fr})
        RadiusThisFrame = mean([Ellipses{fr}(:,3),Ellipses{fr}(:,4)]'); %the 3rd and 4th column contains size info.
        meanFrameRadius = mean(RadiusThisFrame(RadiusThisFrame>0));
    end
        RadiusPerFrame = [RadiusPerFrame meanFrameRadius];
end

% we'll do a dynamic thresholding to prioritize first the schnitz that are very
% close and then incrementally search further away.
Thresholds = 1:0.5:2; %this is a multiplicative factor, 'how many radii' away from the center.
nThresh = size(Thresholds,2);

%initialize information about the stitching state of each schnitz
for j=1:length(schnitzcells)
    schnitzcells(j).ExtendedIntoFutureWithThisThresh = false;
    schnitzcells(j).AlreadyUsed = false; %tells you wether this one has been appended to another one already
    schnitzcells(j).StitchedTo = []; %to know to which schnitz this one was pasted to
    schnitzcells(j).StitchedFrom = []; %to know what original schnitz were pasted to this one to form the final one
end

Original_schnitzcells = schnitzcells; % save the original state before running this script
save([DropboxFolder,filesep,Prefix '_PreStitched.mat'],'Original_schnitzcells');
%% Start stitching loop 

disp('start stitching')
ThingsHaveChangedFlag = true;
IterCounter = 0;
while ThingsHaveChangedFlag
    SomethingToStitchFound = [];
    for thresh = Thresholds %loop over expanding distance thresholds 
        disp(['trying a distance threshold of ' num2str(thresh) 'radii'])
        % reset for a new distance threshold       
        for j=1:length(schnitzcells)
            schnitzcells(j).ExtendedIntoFutureWithThisThresh = false;
        end

        for s1 = 1:length(schnitzcells) 

            if ~(schnitzcells(s1).AlreadyUsed)

                %get info from the last frame of this schnitz
                LastFrame1 = schnitzcells(s1).frames(end);
                FinalXYPos1 = [schnitzcells(s1).cenx(end),schnitzcells(s1).ceny(end)];
                RadiusThisFrame = RadiusPerFrame(LastFrame1); 

                %find schnitz that start the frame after this one ends
                SchnitzThatStartNextFrame = [];
                for s2 = 1:length(schnitzcells)
                    frames = schnitzcells(s2).frames;
                    if frames(1) == LastFrame1+1
                        SchnitzThatStartNextFrame = [SchnitzThatStartNextFrame s2];
                    end
                end

                %find the ones that are also in roughly the same position
                StartNextFrameAndClose = [];
                Distances = [];
                for s3 = SchnitzThatStartNextFrame
                    XYPos = [schnitzcells(s3).cenx(1) schnitzcells(s3).ceny(1)];
                    Distance = pdist([double(FinalXYPos1);double(XYPos)],'euclidean');
                    if  Distance < RadiusThisFrame*(thresh)
                        StartNextFrameAndClose = [StartNextFrameAndClose s3];
                        Distances = [Distances Distance];
                    end
                end

                %if two were close, pick the closest one
                if ~isempty(StartNextFrameAndClose)
                    ClosestSchnitzNextFrame = StartNextFrameAndClose(find(Distances==min(Distances)));
                    ClosestSchnitzNextFrame = ClosestSchnitzNextFrame(1);
                end

                % append the schnitz from the future into the one from the past if
                % we haven't already done it before
                if schnitzcells(s1).ExtendedIntoFutureWithThisThresh == false &&...
                        ~isempty(StartNextFrameAndClose)
                    schnitzcells = AppendSchnitzcellsData(schnitzcells,s1,ClosestSchnitzNextFrame);
                    % mark the schnitz from the past so that nothing else is
                    % appended to it anymore
                    schnitzcells(s1).ExtendedIntoFutureWithThisThresh = true;
                    schnitzcells(s2).AlreadyUsed = true;
                    SomethingToStitchFound = [SomethingToStitchFound 1];
                    length(SomethingToStitchFound);
                end
            end
        end
    end
    ThingsHaveChangedFlag = ~isempty(SomethingToStitchFound);
    IterCounter = IterCounter+1;
    disp(['done with iteration #' num2str(IterCounter) 'of stitching']);
    if ~ThingsHaveChangedFlag
        disp('stitching not changing anymore, exiting loop')
    end
end
disp('done stitching!')


%%
postStitching_schnitzcells = schnitzcells;
save([DropboxFolder,filesep,Prefix '_PostStitch.mat'],'postStitching_schnitzcells');
SchnitzToKillIndices = [schnitzcells.AlreadyUsed];
schnitzcells(SchnitzToKillIndices) = [];
save([DropboxFolder,filesep,Prefix '_lin.mat'],'schnitzcells');


%% Save everything and break schnitzs at mitosis
Stitched_before_breakup = schnitzcells;
save([DropboxFolder,filesep,Prefix '_PreBroken.mat'],'Stitched_before_breakup');
[schnitzcells, Ellipses] = breakUpSchnitzesAtMitoses(schnitzcells, Ellipses, ncVector, nFrames);
[Ellipses, schnitzcells] = addSchnitzIndexToEllipses(Ellipses, schnitzcells);
save2([DropboxFolder,filesep,Prefix '_lin.mat'],schnitzcells);
save2([DropboxFolder,filesep,'Ellipses.mat'],Ellipses);
% TrackNuclei(Prefix,'nWorkers', nWorkers, 'noStitch', 'retrack', 'integrate');

end

