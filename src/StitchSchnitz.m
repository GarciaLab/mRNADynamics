function schnitzcells = StitchSchnitz(varargin)

%This function joins schnitzcells that overlap in space and are contiguous in time.



%Information about about folders
Prefix = varargin{1};

[SourcePath, FISHPath, DefaultDropboxFolder, DropboxFolder, MS2CodePath, PreProcPath,...
configValues, movieDatabasePath] = DetermineAllLocalFolders(Prefix);

FilePrefix=[Prefix,'_'];


[Date, ExperimentType, ExperimentAxis, CoatProtein, StemLoop, APResolution,...
Channel1, Channel2, Objective, Power, DataFolder, DropboxFolderName, Comments,...
nc9, nc10, nc11, nc12, nc13, nc14, CF] = getExperimentDataFromMovieDatabase(Prefix, DefaultDropboxFolder)

ncVector=[0,0,0,0,0,0,0,0,nc9,nc10,nc11,nc12,nc13,nc14];

% Load all the nuclear segmentation and tracking information, and FrameInfo
load([DropboxFolder,filesep,Prefix,filesep,[Prefix '_lin.mat']])
load([DropboxFolder,filesep,Prefix,filesep,'Ellipses.mat'])
load([DropboxFolder,filesep,Prefix,filesep,'FrameInfo.mat'])


%% /\/\/\/\/\/\/*\/*\/*\ END OF LOADING STEPS /*\/*\/*\/\/\/\/\/\/\/\/\/\

%% /|\/|\/|\/|\ Setup stuff /|\/|\/|\/|\

%Start the stitching
[nFrames,~] = size(Ellipses); %how many frames do we have?
Radii = []; %a vector of length = frames that will contain ellipse radius per frame
%get the size information of ellipses in time
%we want to use size info as a threshold to decide if two schnitz in two contiguous frames
%correspond to the same nucleus.
for fr = 1:nFrames
    Radius = Ellipses{fr}(1,3); %the third column contains size info. by definition all ellipses/frame are equal in size
    %Is this value the actual radius? or is it a diameter?
    Radii = [Radii Radius];
end

Thresholds = [1.05:0.05:1.75]; %this number indicates how many radii of distance...
%make two time-contiguous ellipses the same one


%set stitching information
for j=1:length(schnitzcells)
    schnitzcells(j).Valid = 1; %assume all are OK
    schnitzcells(j).StitchedTo = []; %to know to which schnitz this one was pasted to
    schnitzcells(j).StitchedFrom = []; %to know what original schnitz were pasted to this one to form the final one
end

%% /|\/|\/|\/|\/|\ Start stitching loop /|\/|\/|\/|\/|\/|\  
h=waitbar(0,'Stitching broken schnitzs');

for i=1:length(Thresholds)
    
    waitbar(i/length(Thresholds),h);
    display (['Trying threshold', num2str(Thresholds(i))]);
    %length(schnitzcells) for debugging
    threshold = Thresholds(i);
    Rule = 1;
  
    while Rule
        for s1 = 1:length(schnitzcells) %outer loop. 
            %get schnitz basic time/space info
            Frames1 = schnitzcells(s1).frames;
            XPos1 = schnitzcells(s1).cenx;
            YPos1 = schnitzcells(s1).ceny;

            %get the numbers we care about: the extremes.
            FirstFrame1 = Frames1(1);
            LastFrame1 = Frames1(end);
            InitialPosition1 = [XPos1(1),YPos1(1)] ;
            FinalPosition1 = [XPos1(end),YPos1(end)];
            Radius = Radii(LastFrame1);

            for s2 = 1:length(schnitzcells) %inner loop to make all pairwise comparisons
                %[s1, s2]
                %get schnitz basic time/space info
                Frames2 = schnitzcells(s2).frames;
                XPos2 = schnitzcells(s2).cenx;
                YPos2 = schnitzcells(s2).ceny;

                %get the numbers we care about: the extremes.
                FirstFrame2 = Frames2(1);
                LastFrame2 = Frames2(end);

                    InitialPosition2 = [XPos2(1),YPos2(1)] ;
                    FinalPosition2 = [XPos2(end),YPos2(end)];

                %compare schnitz s1 and s2
                if schnitzcells(s1).Valid == 1
                    if LastFrame1+1 == FirstFrame2 && ...
                            pdist([FinalPosition1;InitialPosition2],'euclidean') < Radius*(threshold)
                        % display('overlapping schnitz found')
                        schnitzcells(s2).Valid = 0; %mark so we don't consider it in the future
                        schnitzcells(s2).StitchedTo = [schnitzcells(s2).StitchedTo s1];

                        %paste info
                        schnitzcells(s1).StitchedFrom = [schnitzcells(s1).StitchedFrom s2];
                        schnitzcells(s1).frames = [Frames1 ; Frames2]; %stitch frames
                        schnitzcells(s1).cenx = [schnitzcells(s1).cenx schnitzcells(s2).cenx]; %
                        schnitzcells(s1).ceny = [schnitzcells(s1).ceny schnitzcells(s2).ceny]; %
                        schnitzcells(s1).len = [schnitzcells(s1).len schnitzcells(s2).len]; %
                        schnitzcells(s1).cellno = [schnitzcells(s1).cellno schnitzcells(s2).cellno]; %
                        %schnitzcells(s1).P = [schnitzcells(s1).P schnitzcells(s2).P];
                        %schnitzcells(s1).D = [schnitzcells(s1).D schnitzcells(s2).D]; 
                        %schnitzcells(s1).E = [schnitzcells(s1).E schnitzcells(s2).E];

                        %The Fluo field is only present in input-output
                        %function mode
                        if isfield(schnitzcells,'Fluo')
                            schnitzcells(s1).Fluo = [schnitzcells(s1).Fluo;schnitzcells(s2).Fluo]; %  
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
        %impostors here. We need to pay attention to not screw up the
        %cross-referencing of mother and daughter nuclei.
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
            schnitzcells(k).Valid = 1; %assume all are OK
            schnitzcells(k).StitchedTo = []; %to know to which schnitz this one was pasted to
            schnitzcells(k).StitchedFrom = []; %to know what original schnitz were pasted to this one to form the final one
        end
    end 
    %threshold %(for debugging)
end
close(h)


[schnitzcells, Ellipses] = breakUpSchnitzesAtMitoses(schnitzcells, Ellipses, ncVector, nFrames)
save([DropboxFolder,filesep,Prefix,filesep,Prefix '_lin.mat'],'schnitzcells')
save([DropboxFolder,filesep,Prefix,filesep,'Ellipses.mat'],'Ellipses')


%% Accesory code to check nuclear traces
% % Does the number of schnitz make sense?
% TotalSchnitz = length(schnitzcells);
% StableNuclei = [1,28,55,100,169]; %fill this vector with the frame corresponding to stable nuclei in each nc
% %This vector will be used to count number of ellipses in these frames.
% RealTotalEllipses = 0;
% for sn = 1:length(StableNuclei)
%     StableFrame = StableNuclei(sn);
%     RealTotalEllipses = RealTotalEllipses + length(Ellipses{StableFrame,1});
% end
% RealTotalEllipses
