function [schnitzcells, Ellipses] = StitchSchnitzv2(Prefix, nWorkers)

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
Thresholds = 1:0.5:2.5; %this is a multiplicative factor, 'how many radii' away from the center.
nThresh = size(Thresholds,2);

%initialize information about the stitching state of each schnitz
for j=1:length(schnitzcells)
    schnitzcells(j).AlreadyUsed = false; %tells you wether this one has been appended to another one already
    schnitzcells(j).StitchedTo = []; %to know to which schnitz this one was pasted to
    schnitzcells(j).StitchedFrom = []; %to know what original schnitz were pasted to this one to form the final one
end

Original_schnitzcells = schnitzcells; % save the original state before running this script
save([DropboxFolder,filesep,Prefix '_PreStitched.mat'],'Original_schnitzcells');
%% Start stitching loop 

for thresh = Thresholds %loop over expanding distance thresholds 
    thresh
    RunOutToStitch_withThisThresh = 1;     
    while RunOutToStitch_withThisThresh
        
        FoundSomethingToStitch = 0;
        for s1 = 1:length(schnitzcells) % outer loop, for the 'stitched to' schnitz
                       
            %get info from the last frame of this schnitz
            LastFrame1 = schnitzcells(s1).frames(end);
            FinalXYPos1 = [schnitzcells(s1).cenx(end),schnitzcells(s1).ceny(end)];
            RadiusThisFrame = RadiusPerFrame(LastFrame1); 
           
            % now we'll do all possible pairwise comparisons to the rest of
            % the schnitzs
            for s2 = 1:length(schnitzcells)
                
                FirstFrame2 = schnitzcells(s2).frames(1);          
                InitialXYPos2 = [schnitzcells(s2).cenx(1),schnitzcells(s2).ceny(1)] ;
               
                %compare schnitz s1 and s2 if s2 hasn't been used already
                if schnitzcells(s2).AlreadyUsed == false && LastFrame1+1 == FirstFrame2 && ...
                     pdist([double(FinalXYPos1);double(InitialXYPos2)],'euclidean') < RadiusThisFrame*(thresh)
                    disp('found something to stitch!!')
                    FoundSomethingToStitch = FoundSomethingToStitch+1;
                    schnitzcells(s2).AlreadyUsed = true;                    
                    schnitzcells(s2).StitchedTo = [schnitzcells(s2).StitchedTo s1];
                    schnitzcells(s1).StitchedFrom = [schnitzcells(s1).StitchedFrom s2];
                    %pass fields' info by appending them to the right
                    %of the outer loop schnitz ones                       
                    schnitzcells = AppendSchnitzcellsData(schnitzcells,s1,s2);
                       
                end
            end

        end
        
        if FoundSomethingToStitch == 0
            RunOutToStitch_withThisThresh = false;
        end
        RunOutToStitch_withThisThresh
    end %end of while
end


%%
postStitching_schnitzcells = schnitzcells;
save([DropboxFolder,filesep,Prefix '_PostStitch.mat'],'postStitching_schnitzcells');
FinalSchnitzToKillIndices = [schnitzcells.AlreadyUsed];
schnitzcells(FinalSchnitzToKillIndices) = [];
save([DropboxFolder,filesep,Prefix '_lin.mat'],'schnitzcells');


%% Save everything and break schnitzs at mitosis
Stitched_before_breakup = schnitzcells;
save([DropboxFolder,filesep,Prefix '_PreBroken.mat'],'Stitched_before_breakup');
[schnitzcells, Ellipses] = breakUpSchnitzesAtMitoses(schnitzcells, Ellipses, ncVector, nFrames);
[Ellipses, schnitzcells] = addSchnitzIndexToEllipses(Ellipses, schnitzcells);
save([DropboxFolder,filesep,Prefix '_lin.mat'],'schnitzcells');
save([DropboxFolder,filesep,'Ellipses.mat'],'Ellipses');
%TrackNuclei(Prefix,'nWorkers', nWorkers, 'noStitch', 'retrack', 'integrate');

end

