function schnitzcells = StitchSchnitz(Prefix, nWorkers)

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
        Radius = Ellipses{fr}(1,3); %the third column contains size info. by definition all ellipses/frame are equal in size
    else
        Radius = nan;
    end
        Radii = [Radii Radius];
end

% we'll do a dynamic thresholding to prioritize first the schnitz that are very
% close.
Thresholds = 1.05:0.05:1.75; 
%these numbers indicate how many radii of distance make two time-contiguous ellipses be considered the same one
nThresh = size(Thresholds);

%initialize information about the stitching state of each schnitz
for j=1:length(schnitzcells)
    schnitzcells(j).Valid = true; %assume all schnitz are good to start with
    schnitzcells(j).StitchedTo = []; %to know to which schnitz this one was pasted to
    schnitzcells(j).StitchedFrom = []; %to know what original schnitz were pasted to this one to form the final one
end

%% /|\/|\/|\/|\/|\ Start stitching loop /|\/|\/|\/|\/|\/|\
h=waitbar(0,'Stitching broken schnitzs');

for i=1:nThresh %loop over shnrinking distance thresholds
    
    try waitbar(i/length(Thresholds),h); catch; end
    %     display (['Trying threshold', num2str(Thresholds(i))]);
    %length(schnitzcells) for debugging
    threshold = Thresholds(i);
    Rule = 1;
    
    while Rule
        for s1 = 1:length(schnitzcells) % outer loop, for the 'stitched from' schnitz
            %get schnitz basic time/space info
            %the number '1' in these variable names refer to the outer loop, 'stitched from' schnitz
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
            %Radius = Radii(LastFrame1); 
            Radius = mean([Ellipses{LastFrame1,1}(EllipseIndex,3),Ellipses{LastFrame1,1}(EllipseIndex,4)]);
           
            % now we'll do all possible pairwise comparisons to the rest of
            % the schnitz
            for s2 = 1:length(schnitzcells) %inner loop for the 'stitched to' schnitz
                % get schnitz basic time/space info
                % the number '2' in these variable names refer to the inner loop, 'stitched to' schnitz
                Frames2 = schnitzcells(s2).frames;
                XPos2 = schnitzcells(s2).cenx;
                YPos2 = schnitzcells(s2).ceny;
                
                %get the numbers we care about: first frame.
                FirstFrame2 = Frames2(1);          
                InitialPosition2 = [XPos2(1),YPos2(1)] ;
                
                %compare schnitz s1 and s2
                if schnitzcells(s1).Valid %
                    % if they are contiguous in time and closer than
                    % 'thresh' in space:
                    if LastFrame1+1 == FirstFrame2 && ...
                            pdist([double(FinalPosition1);double(InitialPosition2)],'euclidean') < Radius*(threshold)
                        % display('overlapping schnitz found')
                        schnitzcells(s2).Valid = false; %mark so we don't consider it in the future
                        schnitzcells(s2).StitchedTo = [schnitzcells(s2).StitchedTo s1];
                        
                        %paste fields' info
                        schnitzcells(s1).StitchedFrom = [schnitzcells(s1).StitchedFrom s2];
                        schnitzcells(s1).frames = [Frames1 ; Frames2]; %stitch frames
                        schnitzcells(s1).cenx = [schnitzcells(s1).cenx schnitzcells(s2).cenx]; %
                        schnitzcells(s1).ceny = [schnitzcells(s1).ceny schnitzcells(s2).ceny]; %
                        if isfield(schnitzcells, 'len')
                            schnitzcells(s1).len = [schnitzcells(s1).len schnitzcells(s2).len]; %
                        end
                        schnitzcells(s1).cellno = [schnitzcells(s1).cellno schnitzcells(s2).cellno]; %            
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
            end
        end
        
        %kill impostor schnitz forever
        %create new schnitzcells struct using only valid schnitz from former one
        ImpostorSchnitz = find ([schnitzcells.Valid] == 0);
        if isempty(ImpostorSchnitz)
            Rule = 0; %finish the while loop if can't find more. Continue with next threshold
        end
        
        %Create a new schnitzcells structure. We'll get rid of the
        %impostors here.
        schnitzcells2=schnitzcells;
        
        %Go through all schnitzcells, find the references to schnitz that
        %are one index higher than each ImpostorSchnitz and subtract one.

        for i=1:length(ImpostorSchnitz)
            for s=1:length(schnitzcells2)
                %If this schnitz has a mother or daughters with an index
                %over ImpostorSchnitz(i), then subtract one. If the index
                %is ImpostorSchnitz(i), then leave this empty.
                
                %Parent cell
                if schnitzcells2(s).P>ImpostorSchnitz(i)
                    schnitzcells2(s).P=schnitzcells2(s).P-1;
                elseif schnitzcells2(s).P==ImpostorSchnitz(i)
                    schnitzcells2(s).P=[];
                end
                %Daughter cell E
                if schnitzcells2(s).E>ImpostorSchnitz(i)
                    schnitzcells2(s).E=schnitzcells2(s).E-1;
                elseif schnitzcells2(s).E==ImpostorSchnitz(i)
                    schnitzcells2(s).E=[];
                end
                %Daughter cell D
                if schnitzcells2(s).D>ImpostorSchnitz(i)
                    schnitzcells2(s).D=schnitzcells2(s).D-1;
                elseif schnitzcells2(s).D==ImpostorSchnitz(i)
                    schnitzcells2(s).D=[];
                end
            end
            
            %Now that we have re-indexed all schnitzs, delete the
            %ImpostorSchnitz(i) one.
            schnitzcells2=schnitzcells2([1:ImpostorSchnitz(i)-1,...
                ImpostorSchnitz(i)+1:end]);
            %Finally, we also need to shift all ImpostorSchnitzses by one
            ImpostorSchnitz=ImpostorSchnitz-1;
        end
        
        schnitzcells = schnitzcells2; %replace original by schnitzcells2
        clear schnitzcells2
        for k=1:length(schnitzcells)
            schnitzcells(k).Valid = true; %assume all are OK
            schnitzcells(k).StitchedTo = []; %to know to which schnitz this one was pasted to
            schnitzcells(k).StitchedFrom = []; %to know what original schnitz were pasted to this one to form the final one
        end
    end
    %threshold %(for debugging)
end

try close(h); catch; end


[schnitzcells, Ellipses] = breakUpSchnitzesAtMitoses(schnitzcells, Ellipses, ncVector, nFrames);
save([DropboxFolder,filesep,Prefix '_lin.mat'],'schnitzcells', '-append');
save([DropboxFolder,filesep,'Ellipses.mat'],'Ellipses', '-append');
TrackNuclei(Prefix,'nWorkers', nWorkers, 'noStitch', 'retrack', 'integrate');

