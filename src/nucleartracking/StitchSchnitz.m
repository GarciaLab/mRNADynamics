function [schnitzcells, Ellipses] = StitchSchnitz(Prefix, nWorkers)

%This function joins schnitzcells that overlap in space and are contiguous in time.

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

%% /\/\/\/\/\/\/*\/*\/*\ END OF LOADING STEPS /*\/*\/*\/\/\/\/\/\/\/\/\/\

%% /|\/|\/|\/|\ Setup stuff /|\/|\/|\/|\

%Start the stitching
[nFrames,~] = size(Ellipses); %how many frames do we have?
Radii = []; %a vector of length = frames that will contain ellipse radius per frame
%get the size information of ellipses in time
%we want to use size info as a threshold to decide if two schnitz in two contiguous frames
%correspond to the same nucleus.
for fr = 1:nFrames
    if ~isempty(Ellipses{fr})
        Radius = mean(Ellipses{fr}(:,3)); %the third column contains size info. by definition all ellipses/frame are equal in size
        Radius = Radius(Radius>0);
    end
        Radii = [Radii Radius];
end

% we'll do a dynamic thresholding to prioritize first the schnitz that are very
% close.
Thresholds = 1:0.5:2.5; 
%these numbers indicate how many radii of distance make two time-contiguous ellipses be considered the same one
nThresh = size(Thresholds,2);

%initialize information about the stitching state of each schnitz
for j=1:length(schnitzcells)
    schnitzcells(j).Valid = true; %assume all schnitz are good to start with
    schnitzcells(j).StitchedTo = []; %to know to which schnitz this one was pasted to
    schnitzcells(j).StitchedFrom = []; %to know what original schnitz were pasted to this one to form the final one
end

Original_schnitzcells = schnitzcells;
save([DropboxFolder,filesep,Prefix '_PreStitched.mat'],'Original_schnitzcells');
%% /|\/|\/|\/|\/|\ Start stitching loop /|\/|\/|\/|\/|\/|\


h=waitbar(0,'Stitching broken schnitzs');

for i=1:nThresh %loop over expanding distance thresholds
    
    try waitbar(i/length(Thresholds),h); catch; end
    %     display (['Trying threshold', num2str(Thresholds(i))]);
    %length(schnitzcells) for debugging
    threshold = Thresholds(i);
    Rule = 1;
    
    while Rule
        for s1 = 1:length(schnitzcells) % outer loop, for the 'stitched to' schnitz
            %get schnitz basic time/space info
            %the number '1' in these variable names refer to the outer loop, 'stitched to' schnitz
            Frames1 = schnitzcells(s1).frames;
            XPos1 = schnitzcells(s1).cenx;
            YPos1 = schnitzcells(s1).ceny;
            
            %get the numbers we care about from the last frame.
            LastFrame1 = Frames1(end);
            FinalPosition1 = [XPos1(end),YPos1(end)];
            % next, assuming that the radius of this schnitz in this frame is equal to 
            % the mean radii of all ellipses in this frame:
            EllipseIndex = schnitzcells(s1).cellno;
            EllipseIndex = EllipseIndex(end);
            Radius = Radii(LastFrame1); 
            %Radius = mean([Ellipses{LastFrame1}(EllipseIndex,3),Ellipses{LastFrame1}(EllipseIndex,4)]);
           
            % now we'll do all possible pairwise comparisons to the rest of
            % the schnitz
            for s2 = 1:length(schnitzcells) %inner loop for the 'stitched from' schnitz
                % get schnitz basic time/space info
                % the number '2' in these variable names refer to the inner loop, 'stitched from' schnitz
                Frames2 = schnitzcells(s2).frames;
                XPos2 = schnitzcells(s2).cenx;
                YPos2 = schnitzcells(s2).ceny;
                
                %get the numbers we care about: the first frame.
                FirstFrame2 = Frames2(1);          
                InitialPosition2 = [XPos2(1),YPos2(1)] ;
                
                %compare schnitz s1 and s2
                if schnitzcells(s1).Valid %
                    % if they are contiguous in time and closer than
                    % radius scaled by 'threshold' in space:
                    if LastFrame1+1 == FirstFrame2 && ...
                            pdist([double(FinalPosition1);double(InitialPosition2)],'euclidean') < Radius*(threshold)
                        %display(num2str(Radius*(threshold)))
                        schnitzcells(s2).Valid = false; %mark so we don't consider it in the future
                        schnitzcells(s2).StitchedTo = [schnitzcells(s2).StitchedTo s1];
                        
                        %pass fields' info by concatenating them to the
                        %outer loop schnitz ones
                        schnitzcells(s1).StitchedFrom = [schnitzcells(s1).StitchedFrom s2];
                        schnitzcells(s1).frames = [Frames1 ; Frames2]; 
                        schnitzcells(s1).cenx = [schnitzcells(s1).cenx schnitzcells(s2).cenx]; 
                        schnitzcells(s1).ceny = [schnitzcells(s1).ceny schnitzcells(s2).ceny];
                        
                        if isfield(schnitzcells, 'len')
                            schnitzcells(s1).len = [schnitzcells(s1).len schnitzcells(s2).len]; 
                        end
                        
                        schnitzcells(s1).cellno = [schnitzcells(s1).cellno schnitzcells(s2).cellno];
                        
                        if isfield(schnitzcells,'Fluo')
                            schnitzcells(s1).Fluo = [schnitzcells(s1).Fluo;schnitzcells(s2).Fluo]; %
                        end                
                        if isfield(schnitzcells, 'APpos')
                            schnitzcells(s1).APpos = [schnitzcells(s1).APpos;schnitzcells(s2).APpos];
                        end
                        if isfield(schnitzcells, 'APpos')
                            schnitzcells(s1).DVpos =  [schnitzcells(s1).DVpos;schnitzcells(s2).DVpos];
                        end
                        if isfield(schnitzcells, 'FrameApproved')
                            schnitzcells(s1).FrameApproved =  [schnitzcells(s1).FrameApproved,schnitzcells(s2).FrameApproved];
                        end
                        if isfield(schnitzcells, 'Approved')
                            schnitzcells(s1).Approved =  [schnitzcells(s1).Approved || schnitzcells(s2).Approved];
                        end
                        if isfield(schnitzcells, 'FluoTimeTrace')
                            schnitzcells(s1).FluoTimeTrace =  [schnitzcells(s1).FluoTimeTrace;schnitzcells(s2).FluoTimeTrace];    
                        end                       
                    end
                end
            end  %end of inner loop
        end  %end of outer loop
        % by now we have 'valid' schnitz, which are the ones that have gotten extended into the future by being concatenated to others 
        % and also 'not valid' schnitz, the ones whose info has been added to a valid one. These we must get rid of:
        ImpostorSchnitz = find ([schnitzcells.Valid] == 0);
        %If we couldn't find any shcnitz to stitch, finish the while loop and continue with the next threshold
        if isempty(ImpostorSchnitz)
            Rule = 0;
        end
        % But if there were stitched schnitz we have to permanently delete them.
        % Create a new schnitzcells structure. We'll get rid of the
        % impostors here.
        schnitzcells2=schnitzcells;
        
        %Go through all schnitzcells, find the references to schnitz that
        %are one index higher than each ImpostorSchnitz and subtract one.

        for k=1:length(ImpostorSchnitz)
            %Now that we have re-indexed all schnitzs, delete the
            %ImpostorSchnitz(i) one. We do this by selecting all the rows
            %before and after the 'ith' one, without including it.
            schnitzcells2 = schnitzcells2([1:ImpostorSchnitz(k)-1,ImpostorSchnitz(k)+1:end]);          
            % We also need to shift all ImpostorSchnitzses by one since in
            % the new schnitzcells all indices after this impostor schnitz
            % got reduced by 1 after we deleted it.
            ImpostorSchnitz=ImpostorSchnitz-1;
        end
        
        schnitzcells = schnitzcells2; %replace original by schnitzcells2
        clear schnitzcells2
        % re-initialize the schnitzcell struct to start all over again
        for k=1:length(schnitzcells)
            schnitzcells(k).Valid = true; %assume all are OK again before using the next, slightly larger distance threshold
%            schnitzcells(k).StitchedTo = []; %to know to which schnitz this one was pasted to
%            schnitzcells(k).StitchedFrom = []; %to know what original schnitz were pasted to this one to form the final one
        end
    end
end %end of thresholds loop

try close(h); catch; end

%% Save everything and break schnitzs at mitosis
Stitched_before_breakup = schnitzcells;
save([DropboxFolder,filesep,Prefix '_PreStitched.mat'],'Stitched_before_breakup');

%[schnitzcells, Ellipses] = breakUpSchnitzesAtMitoses(schnitzcells, Ellipses, ncVector, nFrames);
[Ellipses, schnitzcells] = addSchnitzIndexToEllipses(Ellipses, schnitzcells);
save([DropboxFolder,filesep,Prefix '_lin.mat'],'schnitzcells');
save([DropboxFolder,filesep,'Ellipses.mat'],'Ellipses');
%TrackNuclei(Prefix,'nWorkers', nWorkers, 'noStitch', 'retrack', 'integrate');

end

