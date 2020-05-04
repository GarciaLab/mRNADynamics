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
    schnitzcells(j).AlreadyUsed = false; %tells you wether this one has been appended to another one already
    schnitzcells(j).StitchedTo = []; %to know to which schnitz this one was pasted to
    schnitzcells(j).StitchedFrom = []; %to know what original schnitz were pasted to this one to form the final one
end

Original_schnitzcells = schnitzcells; % save the original state before running this script
save([DropboxFolder,filesep,Prefix '_PreStitched.mat'],'Original_schnitzcells');
%% Start stitching loop 

disp('stitching')

for thresh = Thresholds %loop over expanding distance thresholds 
    h = waitbar(0,['Stitching broken schnitzs using ' num2str(thresh) ' radii of distance']);
    RunOutToStitch_withThisThresh = 1;     
    while RunOutToStitch_withThisThresh
        
        FoundSomethingToStitch = 0;
        for s1 = 1:length(schnitzcells) % outer loop, for the 'stitched to' schnitz
            try waitbar(s1/length(schnitzcells),h); catch; end
           
            %get info from the last frame of this schnitz
            LastFrame1 = schnitzcells(s1).frames(end);
            FinalXYPos1 = [schnitzcells(s1).cenx(end),schnitzcells(s1).ceny(end)];
            RadiusThisFrame = RadiusPerFrame(LastFrame1); 
           
            SchnitzThatStartNextFrame = [];
            for s2 = 1:length(schnitzcells)
                frames = schnitzcells(s2).frames;
                if frames(1) == LastFrame1+1
                    SchnitzThatStartNextFrame = [SchnitzThatStartNextFrame s2];
                end
            end                        
            % now we'll do all possible pairwise comparisons to the rest of
            % the schnitzs
            for s3 = SchnitzThatStartNextFrame
                InitialXYPos2 = [schnitzcells(s3).cenx(1),schnitzcells(s3).ceny(1)] ;              
                %compare schnitz s1 and s3 if s3 hasn't been used already
                if schnitzcells(s3).AlreadyUsed == false && ...
                     pdist([double(FinalXYPos1);double(InitialXYPos2)],'euclidean') < RadiusThisFrame*(thresh)
%                     %update the last frame info of the s1 schnitz
%                     LastFrame1 = schnitzcells(s3).frames(end);
%                     FinalXYPos1 = [schnitzcells(s3).cenx(end),schnitzcells(s3).ceny(end)];
%                     RadiusThisFrame = RadiusPerFrame(LastFrame1); 
                    
                    FoundSomethingToStitch = FoundSomethingToStitch+1;
                    schnitzcells(s3).AlreadyUsed = true;                    
                    schnitzcells(s3).StitchedTo = [schnitzcells(s3).StitchedTo s1];
                    schnitzcells(s1).StitchedFrom = [schnitzcells(s1).StitchedFrom s3];
                    %pass fields' info by appending them to the right
                    %of the outer loop schnitz ones                       
                    schnitzcells = AppendSchnitzcellsData(schnitzcells,s1,s3);                      
                end
            end

        end
        
        if FoundSomethingToStitch == 0
            RunOutToStitch_withThisThresh = false;
        end
    end %end of while
    try close(h); catch; end
end

disp('done stitching')


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
save2([DropboxFolder,filesep,Prefix '_lin.mat'],schnitzcells);
save2([DropboxFolder,filesep,'Ellipses.mat'],Ellipses);
% TrackNuclei(Prefix,'nWorkers', nWorkers, 'noStitch', 'retrack', 'integrate');

end

