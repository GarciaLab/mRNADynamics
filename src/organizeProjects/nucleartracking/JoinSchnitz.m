function [Particles,schnitzcells]=JoinSchnitz(Particles,schnitzcells,TargetSchnitz,SourceSchnitz,CurrentFrame)

%Joins two elements of schnitzcells. We need to be careful about
%renumbering of the nuclear indices in Particles

%V2: Updated to support Laurent's schnitzcelss

%First, check that there isn't another particle associated with one of the
%schnitz
for i=1:length(Particles)
    if ~isempty(Particles(i).Nucleus)
        AssignedNuclei(i)=Particles(i).Nucleus;
    else
        AssignedNuclei(i)=nan;
    end
end


if sum(AssignedNuclei==TargetSchnitz)>1
    error('Problem with the Target schnitz. Is it associated with more than one particle?')
elseif sum(AssignedNuclei==SourceSchnitz)~=0
    %If the source schnitz (the one we clicked on) already has a particle
    %we need to see if we can join both of them.
    
    TargetParticle=find(AssignedNuclei==TargetSchnitz);
    SourceParticle=find(AssignedNuclei==SourceSchnitz);
    
    if ~sum(ismember(Particles(TargetParticle).Frame,Particles(SourceParticle).Frame))
        Particles=JoinParticleTraces(TargetParticle,SourceParticle,Particles);
    else
        %If there is an overlap we will have to split the SourceSchnitz. We
        %do this starting from the frame after CurrentFrame
        [Particles,schnitzcells]=SplitSchnitz(Particles,schnitzcells,...
            CurrentFrame+1,SourceParticle);
        
        %Recalculate the AssignedNuclei information
        for i=1:length(Particles)
            if ~isempty(Particles(i).Nucleus)
                AssignedNuclei(i)=Particles(i).Nucleus;
            else
                AssignedNuclei(i)=nan;
            end
        end
        
        
        
        %error('There is an overlap of frames between the target and source particles. Deal with this case')
    end
        
    

    

    
    %error('Problem with the Source schnitz. is it associated with a particle?')
% else
%     %Find the particles associated with one of the schnitz
%     for i=1:length(Particles)
%         if Particles(i).Nucleus==TargetSchnitz
%             TargetParticle=i;
%         end
%     end

end
   


%Check if there is an overlap of frames. This would mean that we need to
%"disconnect" the source schnitz and potentially the target schnitz
if ~isempty(intersect(schnitzcells(SourceSchnitz).frames,schnitzcells(TargetSchnitz).frames))
    
    %See if we need to connect backwards
    if sum((schnitzcells(TargetSchnitz).frames)<CurrentFrame)==0
        FrameFilterSource=((schnitzcells(SourceSchnitz).frames)<=CurrentFrame);
        FrameFilterTarget=logical(ones(size(schnitzcells(TargetSchnitz).frames)));
    else
        FrameFilterSource=(schnitzcells(SourceSchnitz).frames>=CurrentFrame);
        FrameFilterTarget=(schnitzcells(TargetSchnitz).frames<CurrentFrame);
    end
    

else
    FrameFilterSource=logical(ones(size(schnitzcells(SourceSchnitz).frames)));
    FrameFilterTarget=logical(ones(size(schnitzcells(TargetSchnitz).frames)));
end


%Find out the first and last frame of each schnitz
SourceFrames=[min(schnitzcells(SourceSchnitz).frames(FrameFilterSource)),max(schnitzcells(SourceSchnitz).frames(FrameFilterSource))];
TargetFrames=[min(schnitzcells(TargetSchnitz).frames(FrameFilterTarget)),max(schnitzcells(TargetSchnitz).frames(FrameFilterTarget))];



%Join the schnitzcells entries
Fields=fieldnames(schnitzcells);

%This is the index I'll use in case I need to create a new schnitz as a
%repository
Nschnitz=length(schnitzcells);


for i=1:length(Fields)
    if strcmp(Fields{i},'P')|strcmp(Fields{i},'E')|strcmp(Fields{i},'D')|strcmp(Fields{i},'Approved')
%         if getfield(schnitzcells(TargetSchnitz),Fields{i})==0
%             schnitzcells(TargetSchnitz)=setfield(schnitzcells(TargetSchnitz),Fields{i},...
%                 getfield(schnitzcells(SourceSchnitz),Fields{i}));
%         end

    else
        InfoSourceAll=getfield(schnitzcells(SourceSchnitz),Fields{i});
        
%         %I had to do this for the cases where len and ang are not the right
%         %length!
%         if (strcmp(Fields{i},'ang')|strcmp(Fields{i},'len'))&(length(InfoSourceAll)~=length(FrameFilterSource))
%             InfoSourceToExport=InfoSourceAll;
%             InfoSourceToKeep=[]
%         else
            InfoSourceToExport=InfoSourceAll(FrameFilterSource);
            InfoSourceToKeep=InfoSourceAll(~FrameFilterSource);
        %end
      
        InfoTargetAll=getfield(schnitzcells(TargetSchnitz),Fields{i});
        InfoTargetToKeep=InfoTargetAll(FrameFilterTarget);
        InfoTargetToMove=InfoTargetAll(~FrameFilterTarget);

        %Delete the information to move from the source schnitz
        schnitzcells(SourceSchnitz)=setfield(schnitzcells(SourceSchnitz),Fields{i},...
            InfoSourceToKeep);
        
        %If present, move the information out of the target schnitz
        if sum(~FrameFilterTarget)
            schnitzcells(Nschnitz+1)=setfield(schnitzcells(TargetSchnitz),Fields{i},...
                InfoTargetToMove);
            schnitzcells(TargetSchnitz)=setfield(schnitzcells(TargetSchnitz),Fields{i},...
                InfoTargetToKeep);
        end
           
        %Move the information from the source schnitz to the target
        %schnitz. I need to check wether I have a row or column vector
        if (size(getfield(schnitzcells(TargetSchnitz),Fields{i}),1)==1)&...
                (size(InfoSourceToExport,1)==1)
            schnitzcells(TargetSchnitz)=setfield(schnitzcells(TargetSchnitz),Fields{i},...
                [getfield(schnitzcells(TargetSchnitz),Fields{i}),InfoSourceToExport]);        
        elseif (size(getfield(schnitzcells(TargetSchnitz),Fields{i}),2)==1)&...
                (size(InfoSourceToExport,2)==1)
            schnitzcells(TargetSchnitz)=setfield(schnitzcells(TargetSchnitz),Fields{i},...
                [getfield(schnitzcells(TargetSchnitz),Fields{i});InfoSourceToExport]);   
        elseif (size(getfield(schnitzcells(TargetSchnitz),Fields{i}),1)==1)&...
                (size(InfoSourceToExport,1)>1)  
            schnitzcells(TargetSchnitz)=setfield(schnitzcells(TargetSchnitz),Fields{i},...
                [getfield(schnitzcells(TargetSchnitz),Fields{i});InfoSourceToExport]); 
        elseif (size(getfield(schnitzcells(TargetSchnitz),Fields{i}),1)>1)&...
                (size(InfoSourceToExport,1)==1)     
            schnitzcells(TargetSchnitz)=setfield(schnitzcells(TargetSchnitz),Fields{i},...
                [getfield(schnitzcells(TargetSchnitz),Fields{i});InfoSourceToExport]); 
            
        end
            
        
    end
end




%Sort the new schnitz according to frames
[schnitzcells(TargetSchnitz).frames,Permutations]=sort(schnitzcells(TargetSchnitz).frames);
schnitzcells(TargetSchnitz).cenx=schnitzcells(TargetSchnitz).cenx(Permutations);
schnitzcells(TargetSchnitz).ceny=schnitzcells(TargetSchnitz).ceny(Permutations);
schnitzcells(TargetSchnitz).cellno=schnitzcells(TargetSchnitz).cellno(Permutations);

%I had to add this because of inconsistencies somewhere in the schnitzcells
%structure
% if length(schnitzcells(TargetSchnitz).ang)==length(Permutations)
%     schnitzcells(TargetSchnitz).ang=schnitzcells(TargetSchnitz).ang(Permutations);
% end

if length(schnitzcells(TargetSchnitz).len)==length(Permutations)
    schnitzcells(TargetSchnitz).len=schnitzcells(TargetSchnitz).len(Permutations);

end


% 
% %Check that the resulting schnitz is not self-referential
% if schnitzcells(TargetSchnitz).D==SourceSchnitz
%     schnitzcells(TargetSchnitz).D=0;
% end
% if schnitzcells(TargetSchnitz).E==SourceSchnitz
%     schnitzcells(TargetSchnitz).E=0;
% end
% 
% 
% %Change the mother and daughter cell information accordingly
% if schnitzcells(SourceSchnitz).P~=0
%     if schnitzcells(schnitzcells(SourceSchnitz).P).D==SourceSchnitz
%         schnitzcells(schnitzcells(SourceSchnitz).P).D=TargetSchnitz;
%     elseif schnitzcells(schnitzcells(SourceSchnitz).P).E==SourceSchnitz
%         schnitzcells(schnitzcells(SourceSchnitz).P).E=TargetSchnitz;
%     end
% end
% 
% 
% 
% 
% if schnitzcells(SourceSchnitz).D~=0
%      schnitzcells(schnitzcells(SourceSchnitz).D).P=TargetSchnitz;
% end
% 
% if schnitzcells(SourceSchnitz).E~=0
%      schnitzcells(schnitzcells(SourceSchnitz).E).P=TargetSchnitz;
% end
% 

%Depending on which one of the two came first copy the corresponding mother
%or daughter nuclei information
if TargetFrames(2)<SourceFrames(1)
    schnitzcells(TargetSchnitz).E=schnitzcells(SourceSchnitz).E;
    schnitzcells(TargetSchnitz).D=schnitzcells(SourceSchnitz).D;
elseif SourceFrames(2)<TargetFrames(1)
    schnitzcells(TargetSchnitz).P=schnitzcells(SourceSchnitz).P;
else
    warning('Am I forgetting a case here?')
end
    

%Make sure this schnitz as not approved
schnitzcells(TargetSchnitz).Approved=0;


%Now, we need to get rid of that schnitz. If we got rid of it completely.
%This is not trivial because
%we'd have to go through all the different Particles and schnitz and
%adjust the structure accordingly. An alternative approach, which we'll
%try here, is just to empty the Source schnitz, as nobody will be
%referring to it.

%Do I need to do this? Isn't this done in the previous lines?
if ~sum(FrameFilterSource)
    Fields=fieldnames(schnitzcells);
    for i=1:length(Fields)
        schnitzcells(SourceSchnitz)=setfield(schnitzcells(SourceSchnitz),Fields{i},...
            []);
    end   
end


    
%end







