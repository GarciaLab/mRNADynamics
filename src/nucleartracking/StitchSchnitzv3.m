function [schnitzcells, Ellipses] = StitchSchnitzv3(Prefix, nWorkers)

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
for frame = 1:nFrames
    if ~isempty(Ellipses{frame})
        RadiusThisFrame = mean([Ellipses{frame}(:,3),Ellipses{frame}(:,4)]'); %the 3rd and 4th column contains size info.
        meanFrameRadius = mean(RadiusThisFrame(RadiusThisFrame>0));
    end
    RadiusPerFrame = [RadiusPerFrame meanFrameRadius];
end

%initialize information about the stitching state of each schnitz
for j=1:length(schnitzcells)
    schnitzcells(j).AlreadyUsed = false;
    schnitzcells(j).ExtendedIntoFutureAlready = false;
    schnitzcells(j).StitchedTo = []; %to know to which schnitz this one was pasted to
    schnitzcells(j).StitchedFrom = []; %to know what original schnitz were pasted to this one to form the final one
    schnitzcells(j).cenx = double(schnitzcells(j).cenx);
    schnitzcells(j).ceny = double(schnitzcells(j).ceny);
end

Original_schnitzcells = schnitzcells; % save the original state before running this script
save([DropboxFolder,filesep,Prefix '_PreStitched.mat'],'Original_schnitzcells');

nSchnitz = length(schnitzcells);
%% Make sure the originial schintzcells struct is good:
% check that all the time-dependent fields have the same length as the
% number of frames
for OGSchnitz = 1:length(schnitzcells)
    FrameLength = numel([schnitzcells(OGSchnitz).frames]);
    assert(length(schnitzcells(OGSchnitz).cellno) == FrameLength);
    assert(length(schnitzcells(OGSchnitz).cenx) == FrameLength);
    assert(length(schnitzcells(OGSchnitz).ceny) == FrameLength);
    assert(length(schnitzcells(OGSchnitz).len) == FrameLength);
end
% disp('schnitzcells struct is good entering stitching')

%% Start stitching loop
thresh = 2; %radii units
disp('Stitching schnitzes...')
ThingsHaveChangedFlag = true;
IterCounter = 0;
while ThingsHaveChangedFlag
    
    for s0 = 1:nSchnitz
        schnitzcells(s0).ExtendedIntoFutureAlready = false;
    end
    
    SomethingToStitchFound = [];
    
    for s1 = 1:nSchnitz
        
        if ~(schnitzcells(s1).ExtendedIntoFutureAlready)
            
            %get info from the last frame of this schnitz
            LastFrame1 = schnitzcells(s1).frames(end);
            expandedRadiusThisFrame = RadiusPerFrame(LastFrame1) * (thresh);
            %find schnitz that start the frame after this one ends
            SchnitzThatBeginNextFrame = [];
            
            for s2 = 1:nSchnitz
                if schnitzcells(s2).frames(1) == LastFrame1+1
                    SchnitzThatBeginNextFrame = [SchnitzThatBeginNextFrame s2];
                end
            end
            
            %find the schnitzes that are also within the expanded radius
            BeginNextFrameAndNear = [];
            DistancesToNextFrameSchnitz = [];
            for s3 = SchnitzThatBeginNextFrame
                Distance =  sqrt((schnitzcells(s1).cenx(end) - schnitzcells(s3).cenx(1))^2 +...
                    (schnitzcells(s1).ceny(end) - schnitzcells(s3).ceny(1))^2);
                if  Distance < expandedRadiusThisFrame
                    BeginNextFrameAndNear = [BeginNextFrameAndNear s3];
                    DistancesToNextFrameSchnitz = [DistancesToNextFrameSchnitz Distance];
                end
            end
            
            %if two were close, pick the closest one
            if ~isempty(BeginNextFrameAndNear)
                ClosestSchnitzNextFrame = BeginNextFrameAndNear(DistancesToNextFrameSchnitz==min(DistancesToNextFrameSchnitz));
                ClosestSchnitzNextFrame = ClosestSchnitzNextFrame(1);
            end
            
            % append the schnitz from the future into the one from the past
            if ~isempty(BeginNextFrameAndNear) && ~schnitzcells(s1).ExtendedIntoFutureAlready...
                    && ~schnitzcells(s1).AlreadyUsed &&...
                    ~schnitzcells(ClosestSchnitzNextFrame).AlreadyUsed
                
                assert(schnitzcells(s1).AlreadyUsed == false)
                schnitzcells(s1).ExtendedIntoFutureAlready = true;
                schnitzcells(s1).StitchedFrom = [schnitzcells(ClosestSchnitzNextFrame).StitchedFrom s1];
                SomethingToStitchFound = [SomethingToStitchFound 1];
                schnitzcells = AppendSchnitzcellsData(schnitzcells,s1,ClosestSchnitzNextFrame);
                schnitzcells(ClosestSchnitzNextFrame).AlreadyUsed = true;
                
            end
        end
    end
    
    ThingsHaveChangedFlag = ~isempty(SomethingToStitchFound);
    IterCounter = IterCounter+1;
%     disp(['done with iteration #' num2str(IterCounter) 'of stitching']);
    if ~ThingsHaveChangedFlag
%         disp('stitching not changing anymore, exiting loop')
    end
end
disp('done stitching!')

%% Make sure the post-stitching schintzcells struct is good:
% check that all the time-dependent fields have the same length as the
% number of frames
for OGSchnitz = 1:length(schnitzcells)
    FrameLength = numel([schnitzcells(OGSchnitz).frames]);
    assert(length(schnitzcells(OGSchnitz).cellno) == FrameLength);
    assert(length(schnitzcells(OGSchnitz).cenx) == FrameLength);
    assert(length(schnitzcells(OGSchnitz).ceny) == FrameLength);
    assert(length(schnitzcells(OGSchnitz).len) == FrameLength);
end
% disp('length of schnitzcells fields is consistent after stitching')

%%
postStitching_schnitzcells = schnitzcells;
% save([DropboxFolder,filesep,Prefix '_PostStitch.mat'],'postStitching_schnitzcells');
SchnitzToKillIndices = [schnitzcells.AlreadyUsed];
schnitzcells(SchnitzToKillIndices) = [];
% save([DropboxFolder,filesep,Prefix '_lin.mat'],'schnitzcells');


%% Save everything and break schnitzs at mitosis
Stitched_before_breakup = schnitzcells;
% save([DropboxFolder,filesep,Prefix '_PreBroken.mat'],'Stitched_before_breakup');
[Ellipses, schnitzcells] = addSchnitzIndexToEllipses(Ellipses, schnitzcells);
[schnitzcells, Ellipses] = breakUpSchnitzesAtMitoses(schnitzcells, Ellipses, ncVector, nFrames);
save2([DropboxFolder,filesep,Prefix '_lin.mat'],schnitzcells);
save2([DropboxFolder,filesep,'Ellipses.mat'],Ellipses);

end

